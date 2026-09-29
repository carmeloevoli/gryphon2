#!/usr/bin/env python3
"""Run the two H=4 kpc source-variation ensembles, matched to halos_H4.

Injection energies and indices vary separately, after generating the same
source catalogue for each seed. Configurations implement model.tex exactly.
"""

import argparse
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path

from run_halos import ROOT, check_build, run_ensemble

MODELS = ("variations_energy", "variations_index")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seeds", type=int, default=1000)
    parser.add_argument("--threads", type=int, default=6, help="threads per model")
    parser.add_argument("--executable", type=Path, default=ROOT / "build/runAnisotropy")
    parser.add_argument("--outdir", type=Path, default=ROOT / "runs/hebreaks_variations_uncut")
    parser.add_argument("--reference", type=Path,
                        default=ROOT / "runs/hebreaks_halos_uncut/halos_H4")
    args = parser.parse_args()
    if args.seeds < 1 or args.threads < 1:
        parser.error("--seeds and --threads must be positive")
    args.executable = args.executable.resolve()
    args.outdir = args.outdir.resolve()
    reference = json.loads((args.reference / "halos_manifest.json").read_text())
    check_build(args.executable, sorted(args.executable.parent.glob("libgryphon.*")), reference)
    with ThreadPoolExecutor(max_workers=len(MODELS)) as pool:
        list(pool.map(lambda name: run_ensemble(
            ROOT / f"configs/hebreaks/{name}.ini", args.outdir / name,
            args, "variations_manifest.json"), MODELS))


if __name__ == "__main__":
    main()
