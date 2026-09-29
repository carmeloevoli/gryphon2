#!/usr/bin/env python3
"""Reproduce Figure 5 using the current config-driven C++ flux calculator.

Fixed seeds 0--9, with no spectral-feature selection. The two halo models may
run concurrently, but seeds within each model are serial (shared params.ini).
Hashes prevent mixing old results with changed code or configuration.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
import os
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parents[1]
SEEDS = tuple(range(10))


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run_model(halo: int, args: argparse.Namespace) -> None:
    config = ROOT / f"configs/hebreaks/figure5_H{halo}.ini"
    directory = args.outdir / f"figure5_H{halo}"
    directory.mkdir(parents=True, exist_ok=True)
    libraries = sorted(args.executable.parent.glob("libgryphon.*"))
    if not libraries:
        raise RuntimeError("cannot fingerprint libgryphon beside the executable")
    identity = {
        "config_text": config.read_text(),
        "config_sha256": digest(config),
        "executable_sha256": digest(args.executable),
        "library_sha256": {p.name: digest(p) for p in libraries},
        "seeds": list(SEEDS),
        "selection": "first ten seeds, no spectral-feature selection",
    }
    manifest = directory / "figure5_manifest.json"
    if manifest.exists():
        if json.loads(manifest.read_text()) != identity:
            raise RuntimeError(f"{directory}: build/config changed; use a new --outdir")
    elif any(directory.glob("flux_*.txt")):
        raise RuntimeError(f"{directory}: unverified output already exists")
    else:
        manifest.write_text(json.dumps(identity, indent=2) + "\n")
    env = dict(os.environ, GRYPHON_CR_THREADS=str(args.threads))
    for seed in SEEDS:
        flux = directory / f"flux_{seed:06d}.txt"
        complete = directory / f"complete_{seed:06d}.json"
        if complete.exists():
            saved = json.loads(complete.read_text())
            if not flux.exists() or saved != {"flux_sha256": digest(flux)}:
                raise RuntimeError(f"{directory}: completed seed {seed} changed")
            continue
        with (directory / f"run_{seed:06d}.log").open("w") as log:
            subprocess.run(
                [str(args.executable), str(config), "--seed", str(seed),
                 "--outdir", str(args.outdir)],
                env=env, stdout=log, stderr=subprocess.STDOUT, check=True,
            )
        complete.write_text(json.dumps({"flux_sha256": digest(flux)}) + "\n")
        print(f"H={halo} kpc: seed {seed} complete ({seed + 1}/{len(SEEDS)})", flush=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--executable", type=Path, default=ROOT / "build/runAnisotropy")
    parser.add_argument("--outdir", type=Path, default=ROOT / "runs/hebreaks_figure5_uncut")
    parser.add_argument("--threads", type=int, default=4, help="threads per halo model")
    args = parser.parse_args()
    if args.threads < 1:
        parser.error("--threads must be positive")
    args.executable = args.executable.resolve()
    args.outdir = args.outdir.resolve()
    with ThreadPoolExecutor(max_workers=2) as pool:
        list(pool.map(lambda halo: run_model(halo, args), (2, 4)))


if __name__ == "__main__":
    main()
