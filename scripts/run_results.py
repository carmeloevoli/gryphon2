#!/usr/bin/env python3
"""Resumable overnight production of all nine proton-paper results ensembles.

Uses Python's standard library and the existing C++ build. It never builds C++,
launches plots, publishes figures, or edits the manuscript. The first five
ensembles may reuse verified current-model outputs; legacy results are excluded.
"""

from __future__ import annotations

import argparse
from concurrent.futures import FIRST_COMPLETED, ThreadPoolExecutor, wait
from contextlib import contextmanager
from dataclasses import dataclass
from datetime import datetime, timezone
import fcntl
import json
import math
import os
from pathlib import Path
import shutil
import signal
import subprocess
import sys
import threading
import time

from run_halos import ROOT, check_build, digest

NAMES = ("halos_H2", "halos_H4", "halos_H8", "variations_energy", "variations_index",
         "low_rate", "associations_independent", "associations_mixed", "associations_clustered")
COMMON = {
    "pid": "H", "spiralmodel": "Xie2024", "transportmodel": "PureDiffusion",
    "injectionmodel": "GalacticRandom", "emin": 1000., "emax": 1e6, "esize": 48.,
    "hkpc": 4., "discsizepc": 50., "rgkpc": 20., "sunkpc": 8.5,
    "d0h": 0.42, "e0": 1000., "delta": 0.36, "ddelta": -1., "snrateyr": 0.02,
    "maxtimemyr": 100., "injslope": 2.34, "injemax": 0., "efficiency": 0.1,
    "varyenergy": "false", "varyslope": "false", "syntheticassociations": "false",
}


def parameters(text):
    result = {}
    for line in text.splitlines():
        line = line.split("#", 1)[0].strip()
        if not line:
            continue
        key, value = (s.strip() for s in line.split("=", 1))
        if key in result:
            raise ValueError(f"duplicate parameter: {key}")
        try:
            result[key] = float(value)
        except ValueError:
            result[key] = value
    return result


def expected_parameters(name, profile="legacy"):
    expected = dict(COMMON, simname=name)
    if profile == "wide":
        expected.update(emin=100., esize=65., maxtimemyr=200.)
    elif profile != "legacy":
        raise ValueError(f"unknown production profile: {profile}")
    if name.startswith("halos_"):
        expected["hkpc"] = float(name[-1])
    elif name == "variations_energy":
        expected["varyenergy"] = "true"
    elif name == "variations_index":
        expected.update(varyslope="true", injslopesigma=0.30)
    elif name == "low_rate":
        expected.update(snrateyr=0.002, efficiency=1.)
    elif name.startswith("associations_"):
        expected.update(syntheticassociations="true", associationmembers=100.,
                        associationradiuspc=30., associationvelocitykms=3.,
                        associationfraction={"independent": 0., "mixed": 0.5,
                                             "clustered": 1.}[name.split("_")[1]])
    return expected


def validate_parameters(text, expected):
    actual = parameters(text)
    for key, value in expected.items():
        if actual.get(key) != value:
            raise ValueError(f"{expected['simname']}: {key}={actual.get(key)!r}, expected {value!r}")
    # The original halo/variation configs leave the one-year default implicit.
    if actual.get("timestep", 1.) != 1.:
        raise ValueError("the independent source sampler must use one-year time steps")


def atomic_json(path, value):
    temporary = path.with_name(path.name + ".partial")
    temporary.write_text(json.dumps(value, indent=2) + "\n")
    temporary.replace(path)


def atomic_copy(source, destination):
    temporary = destination.with_name(destination.name + ".partial")
    shutil.copy2(source, temporary)
    temporary.replace(destination)


@contextmanager
def exclusive_run(outdir):
    outdir.mkdir(parents=True, exist_ok=True)
    with (outdir / ".runner.lock").open("a+") as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            raise RuntimeError(f"another results runner is already using {outdir}") from None
        try:
            yield lock.fileno()
        finally:
            fcntl.flock(lock, fcntl.LOCK_UN)


@dataclass
class Model:
    name: str
    config: Path
    executable: Path
    directory: Path
    manifest_name: str
    libraries: list
    identity: dict
    reuse: Path | None
    completed: set
    new: int = 0
    seconds: float = 0.
    batches: int = 0
    profile: str = "legacy"

    @property
    def expected(self):
        return expected_parameters(self.name, self.profile)

    @property
    def stems(self):
        return ("flux", "population") if self.name.startswith("associations_") else ("flux",)


def make_model(name, args):
    associated = name.startswith("associations_")
    profile = getattr(args, "profile", "legacy")
    config_dir = ROOT / "configs/hebreaks"
    if profile == "wide":
        config_dir /= "wide"
    elif associated or name == "low_rate":
        config_dir /= "results"
    config = config_dir / f"{name}.ini"
    text = config.read_text()
    validate_parameters(text, expected_parameters(name, profile))
    executable = args.build_dir / ("runAssociations" if associated else "runAnisotropy")
    if not executable.is_file() or not os.access(executable, os.X_OK):
        raise RuntimeError(f"missing executable {executable}; build it before starting production")
    libraries = sorted(args.build_dir.glob("libgryphon.*"))
    if not libraries:
        raise RuntimeError(f"no libgryphon found in {args.build_dir}")
    identity = {"config_text": text, "config_sha256": digest(config),
                "executable_sha256": digest(executable),
                "library_sha256": {p.name: digest(p) for p in libraries}}
    reuse = None
    if associated:
        manifest = "benchmark_manifest.json"
        identity.update(model=name.split("_")[1], delay_model="Zapartas2017-single-star-Eq10")
    else:
        identity.update(selection="consecutive seeds from zero, no spectral-feature selection",
                        matched_catalogues=name != "low_rate")
        if name.startswith("halos_"):
            manifest = "halos_manifest.json"
            reuse = ROOT / "runs/hebreaks_halos_uncut" / name
        elif name.startswith("variations_"):
            manifest = "variations_manifest.json"
            reuse = ROOT / "runs/hebreaks_variations_uncut" / name
        else:
            manifest = "low_rate_manifest.json"
    if args.no_reuse or profile == "wide":
        reuse = None
    return Model(name, config, executable, args.outdir / name, manifest,
                 libraries, identity, reuse, set(), profile=profile)


def validate_identity(model, directory):
    manifest = directory / model.manifest_name
    if manifest.exists():
        if json.loads(manifest.read_text()) != model.identity:
            raise RuntimeError(f"{directory}: build/config mismatch; use a fresh --outdir "
                               "(and --no-reuse if the earlier build is incompatible)")
    elif any(directory.glob("flux_*.txt")) or any(directory.glob("complete_*.json")):
        raise RuntimeError(f"{directory}: output exists without the expected manifest")


def checked_completion(model, directory, seed):
    marker = directory / f"complete_{seed:06d}.json"
    if not marker.exists():
        return False
    files = {stem: directory / f"{stem}_{seed:06d}.txt" for stem in model.stems}
    if not all(p.is_file() for p in files.values()):
        raise RuntimeError(f"{directory}: completed seed {seed} has missing output")
    hashes = {f"{stem}_sha256": digest(p) for stem, p in files.items()}
    if json.loads(marker.read_text()) != hashes:
        raise RuntimeError(f"{directory}: completed seed {seed} failed checksum verification")
    return True


def read_table(path, name, seed):
    metadata, rows = {}, []
    for line in path.read_text().splitlines():
        if line.startswith("#"):
            if "=" in line:
                key, value = line[1:].split("=", 1)
                metadata[key.strip()] = value.strip()
        elif line.strip():
            rows.append([float(v) for v in line.split()])
    if metadata.get("simname") != name or metadata.get("seed") != str(seed):
        raise RuntimeError(f"{path}: incorrect model/seed header")
    if not rows or any(not math.isfinite(v) for row in rows for v in row):
        raise RuntimeError(f"{path}: empty or nonfinite output")
    return rows


def validate_outputs(model, directory, seed):
    rows = read_table(directory / f"flux_{seed:06d}.txt", model.name, seed)
    expected = model.expected
    bins = int(expected["esize"])
    if len(rows) != bins or any(len(row) != 6 for row in rows):
        raise RuntimeError(f"{directory}: expected {bins} flux rows with six columns")
    for i, row in enumerate(rows):
        energy = expected["emin"] * (expected["emax"] / expected["emin"]) ** (i / (bins - 1))
        if not math.isclose(row[0], energy, rel_tol=1e-8) or row[1] <= 0:
            raise RuntimeError(f"{directory}: invalid flux or energy grid at seed {seed}")
    if "population" in model.stems:
        rows = read_table(directory / f"population_{seed:06d}.txt", model.name, seed)
        if len(rows) != 1 or len(rows[0]) != 5 or any(v < 0 for v in rows[0]):
            raise RuntimeError(f"{directory}: invalid population counts")
        parents, field, clustered, retained, count_expected = rows[0]
        if count_expected != expected["snrateyr"] * expected["maxtimemyr"] * 1e6 or retained > field + clustered:
            raise RuntimeError(f"{directory}: inconsistent population normalization")


def prepare(model, target):
    model.directory.mkdir(parents=True, exist_ok=True)
    manifest = model.directory / model.manifest_name
    if not manifest.exists():
        atomic_json(manifest, model.identity)
    snapshot = model.directory / "input.ini"
    if snapshot.exists() and digest(snapshot) != model.identity["config_sha256"]:
        raise RuntimeError(f"{snapshot}: stored input changed")
    if not snapshot.exists():
        atomic_copy(model.config, snapshot)
    if digest(snapshot) != model.identity["config_sha256"]:
        raise RuntimeError(f"{snapshot}: configuration changed during initialization")
    model.config = snapshot  # Freeze input; edits to repository configs cannot affect this run.
    copied = 0
    for seed in range(target):
        if checked_completion(model, model.directory, seed):
            model.completed.add(seed)
        elif model.reuse and checked_completion(model, model.reuse, seed):
            validate_outputs(model, model.reuse, seed)
            for stem in model.stems:
                filename = f"{stem}_{seed:06d}.txt"
                atomic_copy(model.reuse / filename, model.directory / filename)
            log = f"run_{seed:06d}.log"
            if (model.reuse / log).exists():
                atomic_copy(model.reuse / log, model.directory / log)
            # Never hard-link outputs: the original verified ensembles remain independent.
            atomic_copy(model.reuse / f"complete_{seed:06d}.json",
                        model.directory / f"complete_{seed:06d}.json")
            model.completed.add(seed)
            copied += 1
    if copied:
        params = model.reuse / "params.ini"
        validate_parameters(params.read_text(), model.expected)
        atomic_copy(params, model.directory / "params.ini")
    if model.completed:
        validate_parameters((model.directory / "params.ini").read_text(), model.expected)
    print(f"{model.name}: {len(model.completed):,}/{target:,} verified ({copied:,} imported)", flush=True)


def run_seed(model, seed, args):
    check_build(model.executable, model.libraries, model.identity)
    if digest(model.config) != model.identity["config_sha256"]:
        raise RuntimeError(f"{model.config}: frozen configuration changed")
    env = dict(os.environ, GRYPHON_CR_THREADS=str(args.threads))
    # Prefer the library we fingerprint, including for an alternate build dir.
    loader_var = "DYLD_LIBRARY_PATH" if sys.platform == "darwin" else "LD_LIBRARY_PATH"
    env[loader_var] = str(model.executable.parent) + (
        os.pathsep + env[loader_var] if env.get(loader_var) else "")
    log_path = model.directory / f"run_{seed:06d}.log"
    with log_path.open("w") as log:
        subprocess.run([str(model.executable), str(model.config), "--seed", str(seed),
                        "--outdir", str(model.directory.parent)], cwd=ROOT, env=env,
                       stdout=log, stderr=subprocess.STDOUT, check=True,
                       pass_fds=(args.lock_fd,) if hasattr(args, "lock_fd") else ())
    check_build(model.executable, model.libraries, model.identity)
    if digest(model.config) != model.identity["config_sha256"]:
        raise RuntimeError(f"{model.config}: frozen configuration changed during a seed")
    validate_parameters((model.directory / "params.ini").read_text(), model.expected)
    validate_outputs(model, model.directory, seed)
    atomic_json(model.directory / f"complete_{seed:06d}.json",
                {f"{stem}_sha256": digest(model.directory / f"{stem}_{seed:06d}.txt")
                 for stem in model.stems})
    model.completed.add(seed)


def run_batch(model, seeds, args, stop):
    for seed in seeds:
        if stop.is_set():
            break
        started = time.monotonic()
        run_seed(model, seed, args)
        model.seconds += time.monotonic() - started
        model.new += 1
    model.batches += 1


def save_status(models, args, started, state, error=None):
    status = {"updated_utc": datetime.now(timezone.utc).isoformat(), "pid": os.getpid(),
              "state": state, "elapsed_hours": (time.monotonic() - started) / 3600,
              "target_per_model": args.seeds, "jobs": args.jobs, "threads_per_job": args.threads,
              "error": error, "models": {m.name: {
                  "completed": len(m.completed), "remaining": args.seeds - len(m.completed),
                  "new_this_session": m.new,
                  "seconds_per_new_realization": m.seconds / m.new if m.new else None,
                  "directory": str(m.directory)} for m in models}}
    atomic_json(args.outdir / "status.json", status)
    return status


def production(models, args, stop):
    started, last_report = time.monotonic(), 0.
    if shutil.disk_usage(args.outdir).free < 2 * 1024**3:
        raise RuntimeError("less than 2 GiB free; free disk space before production")
    error = None
    active = {}
    try:
        with ThreadPoolExecutor(max_workers=args.jobs) as pool:
            while True:
                now = time.monotonic()
                if args.max_hours and now - started >= args.max_hours * 3600:
                    stop.set()
                busy = {m.name for m in active.values()}
                # Short, round-robin batches: every scenario makes progress overnight.
                available = sorted((m for m in models if m.name not in busy and
                                    len(m.completed) < args.seeds),
                                   key=lambda m: (m.batches, len(m.completed), NAMES.index(m.name)))
                if not stop.is_set():
                    for model in available[:args.jobs - len(active)]:
                        seeds = [i for i in range(args.seeds) if i not in model.completed][:args.batch_size]
                        active[pool.submit(run_batch, model, seeds, args, stop)] = model
                if not active:
                    break
                done, _ = wait(active, timeout=2, return_when=FIRST_COMPLETED)
                for future in done:
                    model = active.pop(future)
                    try:
                        future.result()
                    except Exception as exc:
                        stop.set()
                        error = f"{model.name}: {exc} (see per-seed run logs in {model.directory})"
                if now - last_report >= 30:
                    save_status(models, args, started, "stopping" if stop.is_set() else "running", error)
                    print(" | ".join(f"{m.name} {len(m.completed):,}/{args.seeds:,}" for m in models), flush=True)
                    last_report = now
                if error:
                    # Stop submitting work; all active models finish at most their current seed.
                    stop.set()
    finally:
        finished = all(len(m.completed) == args.seeds for m in models)
        state = "failed" if error else "complete" if finished else "stopped"
        save_status(models, args, started, state, error)
    if error:
        raise RuntimeError(error)
    print(f"{state.upper()}: {sum(len(m.completed) for m in models):,} verified spectra. "
          f"Outputs: {args.outdir}", flush=True)
    return 0 if finished else 130


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profile", choices=("legacy", "wide"), default="legacy",
                        help="wide: 100 GeV--1 PeV, 65 energies, 200 Myr; never imports legacy runs")
    parser.add_argument("--seeds", type=int, default=10000, help="total per model, not additional seeds")
    parser.add_argument("--jobs", type=int, default=min(6, max(1, (os.cpu_count() or 2) // 2)))
    parser.add_argument("--threads", type=int, default=2, help="C++ threads per concurrent model")
    parser.add_argument("--batch-size", type=int, default=25, help="seeds before rotating to another model")
    parser.add_argument("--outdir", type=Path)
    parser.add_argument("--build-dir", type=Path, default=ROOT / "build")
    parser.add_argument("--models", nargs="+", choices=NAMES, default=list(NAMES))
    parser.add_argument("--no-reuse", action="store_true", help="do not import previous current-model runs")
    parser.add_argument("--dry-run", action="store_true", help="read-only preflight; no simulations or output writes")
    parser.add_argument("--status", action="store_true", help="show the saved progress report and exit")
    parser.add_argument("--max-hours", type=float, default=0., help="optional time limit; 0 means run to completion")
    parser.add_argument("--allow-sleep", action="store_true", help="disable the macOS idle-sleep inhibitor")
    args = parser.parse_args(argv)
    if args.outdir is None:
        args.outdir = ROOT / ("runs/hebreaks_wide_10k" if args.profile == "wide" else "runs/hebreaks_results_10k")
    if (min(args.seeds, args.jobs, args.threads, args.batch_size) < 1 or
            not math.isfinite(args.max_hours) or args.max_hours < 0):
        parser.error("counts must be positive and --max-hours nonnegative")
    if len(set(args.models)) != len(args.models):
        parser.error("do not repeat model names")
    args.outdir, args.build_dir = args.outdir.resolve(), args.build_dir.resolve()
    if args.outdir in (ROOT, ROOT.parent, Path.home(), Path("/")):
        parser.error("--outdir must be a dedicated simulation directory")
    if args.status:
        status_file = args.outdir / "status.json"
        if not status_file.exists():
            parser.error(f"no production status yet at {status_file}")
        saved = json.loads(status_file.read_text())
        if saved["state"] in ("running", "stopping"):
            try:
                os.kill(saved["pid"], 0)
            except ProcessLookupError:
                saved["state"] = "interrupted (saved runner PID no longer exists)"
        print(json.dumps(saved, indent=2))
        return 0
    models = [make_model(name, args) for name in args.models]
    # Preflight every model before starting any simulation or copying outputs.
    for model in models:
        validate_identity(model, model.directory)
        if model.reuse and model.reuse.exists():
            validate_identity(model, model.reuse)
        else:
            model.reuse = None
    print(f"Target: {args.seeds:,} realizations in each of {len(models)} scenarios; "
          f"up to {args.jobs} models x {args.threads} threads.\nOutput: {args.outdir}", flush=True)
    if args.jobs * args.threads > (os.cpu_count() or 1):
        print("Warning: requested threads exceed the reported CPU count.", flush=True)
    for model in models:
        existing = sum((model.directory / f"complete_{i:06d}.json").exists() for i in range(args.seeds))
        reusable = sum((model.reuse / f"complete_{i:06d}.json").exists() and
                       not (model.directory / f"complete_{i:06d}.json").exists()
                       for i in range(args.seeds)) if model.reuse else 0
        print(f"  {model.name}: {existing:,} existing, {reusable:,} eligible for verified import, "
              f"{args.seeds - existing - reusable:,} new", flush=True)
    if args.dry_run:
        print("Preflight only; completion checksums will be verified before reuse.")
        return 0
    stop = threading.Event()
    for signum in (signal.SIGINT, signal.SIGTERM):
        signal.signal(signum, lambda *_: stop.set())
    inhibitor = None
    with exclusive_run(args.outdir) as lock_fd:
        # Children retain the lock if the Python parent is forcibly killed;
        # a restart cannot collide with a still-writing orphaned C++ process.
        args.lock_fd = lock_fd
        # Close the preflight/lock race with another runner finishing initialization.
        for model in models:
            validate_identity(model, model.directory)
            prepare(model, args.seeds)
        plan = {"description": ("Wide cutoff-free model; 200 Myr; 65 energies from 100 GeV to 1 PeV"
                                if args.profile == "wide" else
                                "Current cutoff-free paper model; 100 Myr; 48 bins from 1 TeV to 1 PeV"),
                "profile": args.profile,
                "models": {m.name: {"manifest": str(m.directory / m.manifest_name),
                                     "reused_from": str(m.reuse) if m.reuse else None,
                                     "config": m.identity["config_text"]} for m in models},
                "independence": "shared catalogues only across halos and injection variations; "
                                "association controls share parameters, not event catalogues",
                "paper_outputs": "not modified; fit-based curvature subsection needs no new Monte Carlo"}
        atomic_json(args.outdir / "results_plan.json", plan)
        if sys.platform == "darwin" and not args.allow_sleep and shutil.which("caffeinate"):
            inhibitor = subprocess.Popen(["caffeinate", "-i", "-w", str(os.getpid())])
        try:
            return production(models, args, stop)
        finally:
            if inhibitor:
                inhibitor.terminate()
                inhibitor.wait()


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (OSError, ValueError, RuntimeError, subprocess.SubprocessError) as error:
        print(f"ERROR: {error}", file=sys.stderr, flush=True)
        sys.exit(1)
