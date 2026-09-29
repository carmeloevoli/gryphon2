#!/usr/bin/env python3
"""Three seed-0 participation diagnostics; never modify overnight ensembles."""

import argparse
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
import subprocess

from run_halos import ROOT, check_build, digest
from run_results import expected_parameters, validate_parameters

MODELS = ("halos_H4", "variations_energy", "variations_index")
SEED = 0  # Chosen a priori, never selected by the resulting spectra or N_eff.


def run_model(name, args):
    reference = args.reference / name
    config = reference / "input.ini"
    validate_parameters(config.read_text(), expected_parameters(name))
    original = reference / "flux_000000.txt"
    if json.loads((reference / "complete_000000.json").read_text()) != {"flux_sha256": digest(original)}:
        raise ValueError(f"{name}: reference flux failed completion check")
    manifest_name = "halos_manifest.json" if name == "halos_H4" else "variations_manifest.json"
    original_manifest = json.loads((reference / manifest_name).read_text())
    if original_manifest["config_sha256"] != digest(config):
        raise ValueError(f"{name}: frozen input differs from its production manifest")
    libraries = sorted(args.executable.parent.glob("libgryphon.*"))
    if not libraries:
        raise ValueError("cannot fingerprint diagnostic physics library")
    identity = {
        "seed": SEED, "selection": "fixed seed 0, chosen before inspecting any diagnostic",
        "config_text": config.read_text(), "config_sha256": digest(config),
        "executable_sha256": digest(args.executable),
        "library_sha256": {p.name: digest(p) for p in libraries},
        "reference_flux_sha256": digest(original), "reference_directory": str(reference),
        "reference_manifest_sha256": digest(reference / manifest_name),
        "diagnostic_source_sha256": {p: digest(ROOT / p) for p in (
            "tools/diagnostics/inspectEffectiveSources.cpp", "include/gryphon/core/participation.h")},
    }
    directory = args.outdir / name
    directory.mkdir(parents=True, exist_ok=True)
    manifest = directory / "participation_manifest.json"
    output = directory / "participation_000000.txt"
    completed = directory / "complete_000000.json"
    if manifest.exists():
        if json.loads(manifest.read_text()) != identity:
            raise ValueError(f"{name}: changed diagnostic inputs/build; use a new --outdir")
    elif output.exists() or completed.exists():
        raise ValueError(f"{name}: unverified outputs already exist")
    else:
        manifest.write_text(json.dumps(identity, indent=2) + "\n")
    if completed.exists():
        if json.loads(completed.read_text()) != {"participation_sha256": digest(output)}:
            raise ValueError(f"{name}: changed completed output")
        print(f"{name}: verified saved diagnostic", flush=True)
        return
    check_build(args.executable, libraries, identity)
    with (directory / "run_000000.log").open("w") as log:
        subprocess.run([str(args.executable), str(config), "--seed", str(SEED),
                        "--outdir", str(args.outdir)], stdout=log, stderr=subprocess.STDOUT, check=True)
    check_build(args.executable, libraries, identity)
    completed.write_text(json.dumps({"participation_sha256": digest(output)}) + "\n")
    print(f"{name}: completed diagnostic", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--executable", type=Path, default=ROOT / "build-effective-sources/inspectEffectiveSources")
    parser.add_argument("--reference", type=Path, default=ROOT / "runs/hebreaks_results_10k")
    parser.add_argument("--outdir", type=Path, default=ROOT / "runs/hebreaks_effective_sources")
    args = parser.parse_args()
    for name in ("executable", "reference", "outdir"):
        setattr(args, name, getattr(args, name).resolve())
    # Both nested directions are forbidden: no diagnostic params.ini can
    # overwrite a production model, nor may a model be its own reference.
    if args.outdir == args.reference or args.outdir in args.reference.parents or args.reference in args.outdir.parents:
        parser.error("diagnostic and reference roots must be separate")
    with ThreadPoolExecutor(max_workers=3) as pool:
        list(pool.map(lambda name: run_model(name, args), MODELS))


if __name__ == "__main__":
    main()
