"""Guard resumability and provenance in the shared halo/variation runner."""

from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

from run_halos import run_ensemble


class EnsembleRunnerTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.config = self.root / "model.ini"
        self.config.write_text("simname = example\n")
        executable = self.root / "runAnisotropy"
        executable.write_text("test executable identity\n")
        self.library = self.root / "libgryphon.test"
        self.library.write_text("test library identity\n")
        self.args = SimpleNamespace(executable=executable, outdir=self.root / "runs",
                                    seeds=2, threads=1)
        self.directory = self.args.outdir / "example"

    def fake_process(self, command, **kwargs):
        seed = int(command[command.index("--seed") + 1])
        (self.directory / f"flux_{seed:06d}.txt").write_text(f"flux for seed {seed}\n")

    def run_fixture(self, callback=None):
        with patch("run_halos.subprocess.run", side_effect=callback or self.fake_process) as run:
            run_ensemble(self.config, self.directory, self.args, "variations_manifest.json")
        return run.call_count

    def test_resume_and_extend_only_unfinished_seeds(self):
        self.assertEqual(self.run_fixture(), 2)
        self.assertEqual(self.run_fixture(), 0)
        self.args.seeds = 3
        self.assertEqual(self.run_fixture(), 1)
        self.assertTrue((self.directory / "complete_000002.json").exists())

    def test_reject_changed_completed_flux(self):
        self.run_fixture()
        (self.directory / "flux_000000.txt").write_text("changed output\n")
        with self.assertRaisesRegex(RuntimeError, "completed seed 0 changed"):
            self.run_fixture()

    def test_reject_changed_configuration_or_build(self):
        self.run_fixture()
        self.config.write_text("simname = changed\n")
        with self.assertRaisesRegex(RuntimeError, "build/config changed"):
            self.run_fixture()
        self.config.write_text("simname = example\n")
        self.library.write_text("changed library\n")
        with self.assertRaisesRegex(RuntimeError, "build/config changed"):
            self.run_fixture()

    def test_no_completion_marker_if_build_changes_mid_seed(self):
        def changing_process(command, **kwargs):
            self.fake_process(command, **kwargs)
            self.library.write_text("changed during seed\n")
        with self.assertRaisesRegex(RuntimeError, "build changed during"):
            self.run_fixture(changing_process)
        self.assertFalse((self.directory / "complete_000000.json").exists())


if __name__ == "__main__":
    unittest.main()
