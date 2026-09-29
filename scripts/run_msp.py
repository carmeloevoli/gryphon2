#!/usr/bin/env python3
"""Launch native runMSP realizations; this script contains no population or transport solver."""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import hashlib
import json
from pathlib import Path
import shlex
import subprocess

ROOT = Path(__file__).resolve().parents[1]
STEMS = ("population", "flux", "contributions", "horizon", "checks", "summary", "params")


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run_one(binary, config, output, seed):
    directory = output / f"seed_{seed:06d}" / config.stem
    complete = directory / f"summary_{seed:06d}.txt"
    if complete.exists():
        if not all((directory / f"{stem}_{seed:06d}.txt").is_file() for stem in STEMS):
            raise RuntimeError(f"Incomplete output set in {directory}")
        return f"seed {seed}: reused completed native run"
    if directory.exists() and any(directory.iterdir()):
        raise RuntimeError(f"Partial run in {directory}; use a fresh --output (nothing was deleted)")
    base = directory.parent
    base.mkdir(parents=True, exist_ok=True)
    command = [str(binary), str(config), "--seed", str(seed), "--outdir", str(base)]
    with (base / "run.log").open("w") as log:
        subprocess.run(command, cwd=ROOT, stdout=log, stderr=subprocess.STDOUT, check=True)
    if not complete.exists():
        raise RuntimeError(f"No completion summary in {directory}")
    return f"seed {seed}: completed"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=ROOT / "configs/msp/jelly.ini")
    parser.add_argument("--binary", type=Path, default=ROOT / "build/runMSP")
    parser.add_argument("--output", type=Path, default=ROOT / "runs/msp_ensemble")
    parser.add_argument("--realizations", type=int, default=100)
    parser.add_argument("--seed", type=int, default=20260924)
    parser.add_argument("--jobs", type=int, default=4)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    if args.realizations < 1 or args.jobs < 1 or args.seed < 0:
        parser.error("Require positive realizations/jobs and a nonnegative first seed")
    binary, config, output = (p.resolve() for p in (args.binary, args.config, args.output))
    if not binary.is_file() or not config.is_file():
        parser.error("Missing executable or config; build the runMSP target first")
    seeds = list(range(args.seed, args.seed + args.realizations))
    identity = dict(binary=str(binary), binary_sha256=digest(binary),
                    config=str(config), config_sha256=digest(config))
    # The CMake runner links this shared library; hashing only the small
    # executable would miss a changed injection/transport implementation.
    libraries = sorted(set(binary.parent.glob("libgryphon.dylib")) |
                       set(binary.parent.glob("libgryphon.so*")))
    if not libraries:
        parser.error("Cannot find libgryphon beside the runner; use the CMake build directory")
    identity["libraries_sha256"] = {str(p): digest(p) for p in libraries}
    manifest_path = output / "manifest.json"
    if manifest_path.exists():
        previous = json.loads(manifest_path.read_text())
        if previous["identity"] != identity:
            parser.error("Existing ensemble has different code/config; choose a new --output")
        if previous["seeds"][0] != seeds[0] or len(previous["seeds"]) > len(seeds):
            parser.error("Existing ensemble can only be resumed or extended from the same first seed")
    elif output.exists() and any(output.iterdir()):
        parser.error("Nonempty output without a manifest; choose a new --output")
    if args.dry_run:
        print(f"{len(seeds)} native realizations, {args.jobs} workers, output {output}")
        print(shlex.join([str(binary), str(config), "--seed", str(seeds[0]),
                          "--outdir", str(output / f"seed_{seeds[0]:06d}")]))
        return
    output.mkdir(parents=True, exist_ok=True)
    manifest = dict(engine="gryphon2/runMSP", identity=identity, seeds=seeds,
                    model=config.stem, status="running")
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        futures = [pool.submit(run_one, binary, config, output, seed) for seed in seeds]
        for index, future in enumerate(as_completed(futures), 1):
            print(f"[{index}/{len(seeds)}] {future.result()}", flush=True)
    manifest["status"] = "complete"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    main()
