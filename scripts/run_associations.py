#!/usr/bin/env python3
"""Resumable matched ensembles; parallelize models, never writes within one model."""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
import os
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parents[1]
MODELS = ("independent", "mixed", "clustered")


def positive(value: str) -> int:
    result = int(value)
    if result < 1:
        raise argparse.ArgumentTypeError("must be positive")
    return result


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run_model(model: str, args: argparse.Namespace) -> None:
    config = args.config_dir / f"associations_{model}.ini"
    directory = args.outdir / f"associations_{model}"
    directory.mkdir(parents=True, exist_ok=True)
    manifest = directory / "benchmark_manifest.json"
    libraries = sorted(args.executable.parent.glob("libgryphon.*"))
    if not libraries:
        raise RuntimeError("cannot fingerprint libgryphon beside the executable")
    identity = {"config_sha256": digest(config), "executable_sha256": digest(args.executable),
                "library_sha256": {path.name: digest(path) for path in libraries},
                "config_text": config.read_text(), "model": model,
                "delay_model": "Zapartas2017-single-star-Eq10"}
    if manifest.exists():
        if json.loads(manifest.read_text()) != identity:
            raise RuntimeError(f"{directory}: config/binary/library changed; use a new output directory")
    elif any(directory.glob("flux_*.txt")):
        raise RuntimeError(f"{directory}: existing output has no benchmark manifest")
    else:
        manifest.write_text(json.dumps(identity, indent=2) + "\n")
    env = dict(os.environ, GRYPHON_CR_THREADS=str(args.threads))
    for seed in range(args.start_seed, args.start_seed + args.seeds):
        flux = directory / f"flux_{seed:06d}.txt"
        population = directory / f"population_{seed:06d}.txt"
        completed = directory / f"complete_{seed:06d}.json"
        if completed.exists():
            saved = json.loads(completed.read_text())
            if flux.exists() and population.exists() and saved == {
                "flux_sha256": digest(flux), "population_sha256": digest(population)
            }:
                continue
            raise RuntimeError(f"{directory}: completed seed {seed} is missing or changed")
        with (directory / f"run_{seed:06d}.log").open("w") as log:
            subprocess.run([str(args.executable), str(config), "--seed", str(seed),
                            "--outdir", str(args.outdir)], env=env, stdout=log,
                           stderr=subprocess.STDOUT, check=True)
        completed.write_text(json.dumps({"flux_sha256": digest(flux),
                                         "population_sha256": digest(population)}) + "\n")
        if (seed - args.start_seed + 1) % 16 == 0 or seed == args.start_seed + args.seeds - 1:
            print(f"{model}: {seed - args.start_seed + 1}/{args.seeds}", flush=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seeds", type=positive, default=256)
    parser.add_argument("--start-seed", type=int, default=0)
    parser.add_argument("--jobs", type=positive, default=3, help="concurrent models (maximum 3)")
    parser.add_argument("--threads", type=positive, default=3, help="energy-bin threads per model")
    parser.add_argument("--executable", type=Path, default=ROOT / "build/runAssociations")
    parser.add_argument("--config-dir", type=Path, default=ROOT / "configs/hebreaks")
    parser.add_argument("--outdir", type=Path, default=ROOT / "runs/hebreaks_associations")
    args = parser.parse_args()
    if args.start_seed < 0:
        parser.error("--start-seed must be nonnegative")
    for name in ("executable", "config_dir", "outdir"):
        setattr(args, name, getattr(args, name).resolve())
    with ThreadPoolExecutor(max_workers=min(args.jobs, len(MODELS))) as pool:
        list(pool.map(lambda model: run_model(model, args), MODELS))


if __name__ == "__main__":
    main()
