"""Production preflight, resumption, data validation, and stopping tests."""

import json
from pathlib import Path
import subprocess
import tempfile
from types import SimpleNamespace
import threading
import unittest
from unittest.mock import patch

from run_halos import digest
from run_results import (NAMES, atomic_json, checked_completion, exclusive_run,
                         expected_parameters, make_model, prepare, production,
                         run_seed, validate_identity, validate_outputs, validate_parameters)


class ResultsRunnerTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        build = self.root / "build"
        build.mkdir()
        for name in ("runAnisotropy", "runAssociations", "libgryphon.test"):
            path = build / name
            path.write_text(f"fixture {name}\n")
            path.chmod(0o755)
        self.args = SimpleNamespace(build_dir=build, outdir=self.root / "output",
                                    no_reuse=True, seeds=2, jobs=2, threads=1,
                                    batch_size=1, max_hours=0)
        self.model = make_model("halos_H4", self.args)

    def fake_process(self, command, **kwargs):
        config = Path(command[1])
        name = next(line.split("=", 1)[1].strip() for line in config.read_text().splitlines()
                    if line.startswith("simname"))
        seed = int(command[command.index("--seed") + 1])
        directory = Path(command[command.index("--outdir") + 1]) / name
        (directory / "params.ini").write_text(config.read_text())
        header = f"# simname = {name}\n# seed = {seed}\n"
        rows = "\n".join(f"{1e3 * 10**(3*i/47):.8e} 1e-5 0.1 0.1 0 0" for i in range(48))
        (directory / f"flux_{seed:06d}.txt").write_text(header + rows + "\n")
        if name.startswith("associations_"):
            fraction = expected_parameters(name)["associationfraction"]
            (directory / f"population_{seed:06d}.txt").write_text(
                header + f"10000 {2e6 * (1-fraction)} {2e6 * fraction} 1000000 2000000\n")

    def test_all_nine_current_model_configs(self):
        self.assertEqual(len(NAMES), 9)
        for name in NAMES:
            m = make_model(name, self.args)
            validate_parameters(m.config.read_text(), expected_parameters(name))
        low = make_model("low_rate", self.args)
        self.assertFalse(low.identity["matched_catalogues"])

    def test_resume_does_not_repeat_verified_seeds(self):
        prepare(self.model, self.args.seeds)
        with patch("run_results.subprocess.run", side_effect=self.fake_process) as run:
            result = production([self.model], self.args, threading.Event())
            self.assertEqual(result, 0)
            self.assertEqual(run.call_count, 2)
        resumed = make_model("halos_H4", self.args)
        validate_identity(resumed, resumed.directory)
        prepare(resumed, self.args.seeds)
        with patch("run_results.subprocess.run") as run:
            self.assertEqual(production([resumed], self.args, threading.Event()), 0)
            run.assert_not_called()
        self.args.seeds = 3
        with patch("run_results.subprocess.run", side_effect=self.fake_process) as run:
            self.assertEqual(production([resumed], self.args, threading.Event()), 0)
            self.assertEqual(run.call_count, 1)

    def test_changed_build_and_unverified_directory_are_rejected(self):
        prepare(self.model, 2)
        manifest = self.model.directory / self.model.manifest_name
        old = json.loads(manifest.read_text())
        old["executable_sha256"] = "changed"
        atomic_json(manifest, old)
        with self.assertRaisesRegex(RuntimeError, "build/config mismatch"):
            validate_identity(self.model, self.model.directory)
        manifest.unlink()
        (self.model.directory / "flux_000000.txt").write_text("unverified")
        with self.assertRaisesRegex(RuntimeError, "without the expected manifest"):
            validate_identity(self.model, self.model.directory)

    def test_reuse_copies_verified_data_without_hard_links(self):
        prepare(self.model, 2)
        with patch("run_results.subprocess.run", side_effect=self.fake_process):
            run_seed(self.model, 0, self.args)
        original = self.model.directory
        self.args.outdir = self.root / "imported"
        imported = make_model("halos_H4", self.args)
        imported.reuse = original
        validate_identity(imported, original)
        prepare(imported, 2)
        self.assertEqual(imported.completed, {0})
        source, copy = (p / "flux_000000.txt" for p in (original, imported.directory))
        self.assertEqual(digest(source), digest(copy))
        self.assertNotEqual(source.stat().st_ino, copy.stat().st_ino)
        copy.write_text("changed")
        with self.assertRaisesRegex(RuntimeError, "checksum"):
            checked_completion(imported, imported.directory, 0)
        self.assertTrue(checked_completion(self.model, original, 0))

    def test_failed_seed_gets_no_marker(self):
        prepare(self.model, 2)
        with patch("run_results.subprocess.run", side_effect=subprocess.CalledProcessError(1, "fixture")):
            with self.assertRaises(subprocess.CalledProcessError):
                run_seed(self.model, 0, self.args)
        self.assertFalse((self.model.directory / "complete_000000.json").exists())

    def test_build_changes_during_seed_are_rejected(self):
        prepare(self.model, 2)
        def changed(command, **kwargs):
            self.fake_process(command, **kwargs)
            self.model.libraries[0].write_text("new library")
        with patch("run_results.subprocess.run", side_effect=changed):
            with self.assertRaisesRegex(RuntimeError, "build changed during"):
                run_seed(self.model, 0, self.args)
        self.assertFalse((self.model.directory / "complete_000000.json").exists())

    def test_incomplete_flux_is_not_accepted(self):
        prepare(self.model, 2)
        def truncated(command, **kwargs):
            self.fake_process(command, **kwargs)
            (self.model.directory / "flux_000000.txt").write_text(
                "# simname = halos_H4\n# seed = 0\n1000 1 0 0 0 0\n")
        with patch("run_results.subprocess.run", side_effect=truncated):
            with self.assertRaisesRegex(RuntimeError, "48 flux rows"):
                run_seed(self.model, 0, self.args)
        self.assertFalse((self.model.directory / "complete_000000.json").exists())

    def test_association_counts_are_validated_and_checksummed(self):
        model = make_model("associations_clustered", self.args)
        prepare(model, 2)
        with patch("run_results.subprocess.run", side_effect=self.fake_process):
            run_seed(model, 0, self.args)
        self.assertTrue(checked_completion(model, model.directory, 0))
        pop = model.directory / "population_000000.txt"
        pop.write_text(pop.read_text().replace("2000000\n", "600000\n"))
        with self.assertRaisesRegex(RuntimeError, "normalization"):
            validate_outputs(model, model.directory, 0)
        with self.assertRaisesRegex(RuntimeError, "checksum"):
            checked_completion(model, model.directory, 0)

    def test_stop_is_resumable(self):
        prepare(self.model, 2)
        stop = threading.Event()
        stop.set()
        with patch("run_results.subprocess.run") as run:
            self.assertEqual(production([self.model], self.args, stop), 130)
            run.assert_not_called()
        status = json.loads((self.args.outdir / "status.json").read_text())
        self.assertEqual(status["state"], "stopped")
        self.assertEqual(status["models"]["halos_H4"]["completed"], 0)

    def test_failure_stops_other_work_and_reports_error(self):
        prepare(self.model, 2)
        stop = threading.Event()
        with patch("run_results.subprocess.run", side_effect=RuntimeError("fixture failed")):
            with self.assertRaisesRegex(RuntimeError, "fixture failed"):
                production([self.model], self.args, stop)
        self.assertTrue(stop.is_set())
        self.assertEqual(json.loads((self.args.outdir / "status.json").read_text())["state"], "failed")

    def test_second_runner_is_locked_out(self):
        with exclusive_run(self.args.outdir):
            with self.assertRaisesRegex(RuntimeError, "already using"):
                with exclusive_run(self.args.outdir):
                    pass
        with exclusive_run(self.args.outdir):
            pass


if __name__ == "__main__":
    unittest.main()
