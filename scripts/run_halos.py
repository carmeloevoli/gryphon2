#!/usr/bin/env python3
"""Run matched, cutoff-free H=2,4,8 kpc ensembles for the proton paper.

Seeds are unselected and shared across halo models. Each model runs serially
(the C++ driver writes a shared params.ini), with models in parallel. Output
hashes and a build/config manifest make interruptions safely resumable and
prevent mixing different models. Increasing --seeds extends the same sample.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
import os
from pathlib import Path
import subprocess
import time

ROOT = Path(__file__).resolve().parents[1]
HALOS = (2, 4, 8)


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def check_build(executable: Path, libraries: list[Path], identity: dict) -> None:
    if (digest(executable) != identity["executable_sha256"] or
            {p.name: digest(p) for p in libraries} != identity["library_sha256"]):
        raise RuntimeError("build changed during the ensemble; do not rebuild a running calculation")


def run_model(halo: int, args: argparse.Namespace) -> None:
    config = ROOT / f"configs/hebreaks/halos_H{halo}.ini"
    directory = args.outdir / f"halos_H{halo}"
    run_ensemble(config, directory, args, "halos_manifest.json")


def run_ensemble(config: Path, directory: Path, args: argparse.Namespace,
                 manifest_name: str) -> None:
    """Shared resumable runner; never execute seeds concurrently in one model."""
    directory.mkdir(parents=True, exist_ok=True)
    libraries = sorted(args.executable.parent.glob("libgryphon.*"))
    if not libraries:
        raise RuntimeError("cannot fingerprint libgryphon beside the executable")
    identity = {
        "config_text": config.read_text(),
        "config_sha256": digest(config),
        "executable_sha256": digest(args.executable),
        "library_sha256": {p.name: digest(p) for p in libraries},
        "selection": "consecutive seeds from zero, no spectral-feature selection",
        "matched_catalogues": True,
    }
    manifest = directory / manifest_name
    if manifest.exists():
        if json.loads(manifest.read_text()) != identity:
            raise RuntimeError(f"{directory}: build/config changed; use a new --outdir")
    elif any(directory.glob("flux_*.txt")):
        raise RuntimeError(f"{directory}: unverified output already exists")
    else:
        manifest.write_text(json.dumps(identity, indent=2) + "\n")
    env = dict(os.environ, GRYPHON_CR_THREADS=str(args.threads))
    started = time.monotonic()
    new = 0
    for seed in range(args.seeds):
        flux = directory / f"flux_{seed:06d}.txt"
        complete = directory / f"complete_{seed:06d}.json"
        if complete.exists():
            saved = json.loads(complete.read_text())
            if not flux.exists() or saved != {"flux_sha256": digest(flux)}:
                raise RuntimeError(f"{directory}: completed seed {seed} changed")
            continue
        check_build(args.executable, libraries, identity)
        with (directory / f"run_{seed:06d}.log").open("w") as log:
            subprocess.run(
                [str(args.executable), str(config), "--seed", str(seed),
                 "--outdir", str(args.outdir)],
                env=env, stdout=log, stderr=subprocess.STDOUT, check=True,
            )
        check_build(args.executable, libraries, identity)
        complete.write_text(json.dumps({"flux_sha256": digest(flux)}) + "\n")
        new += 1
        if seed == 0 or (seed + 1) % 25 == 0 or seed + 1 == args.seeds:
            elapsed = time.monotonic() - started
            print(f"{directory.name}: {seed + 1}/{args.seeds} complete; "
                  f"{elapsed / new:.2f} s/new realization", flush=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seeds", type=int, default=1000)
    parser.add_argument("--executable", type=Path, default=ROOT / "build/runAnisotropy")
    parser.add_argument("--outdir", type=Path, default=ROOT / "runs/hebreaks_halos_uncut")
    parser.add_argument("--threads", type=int, default=4, help="threads per halo model")
    args = parser.parse_args()
    if args.seeds < 1 or args.threads < 1:
        parser.error("--seeds and --threads must be positive")
    args.executable = args.executable.resolve()
    args.outdir = args.outdir.resolve()
    with ThreadPoolExecutor(max_workers=len(HALOS)) as pool:
        list(pool.map(lambda halo: run_model(halo, args), HALOS))


if __name__ == "__main__":
    main()
