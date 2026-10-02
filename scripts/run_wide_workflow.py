"""Run wide production, then automatically verify and analyse all nine ensembles.

The preview and final outputs live in the run directory. No manuscript files
are overwritten. Interrupted production resumes through run_results.py.
"""
import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import signal
import subprocess
import sys

from run_halos import ROOT, digest
from run_results import atomic_json, exclusive_run


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path, default=ROOT / "runs/hebreaks_wide_10k")
    parser.add_argument("--build-dir", type=Path, default=ROOT / "build-hebreaks-wide")
    parser.add_argument("--history-check", type=Path,
                        default=ROOT / "runs/hebreaks_history_check/history_H8/convergence.json")
    parser.add_argument("--seeds", type=int, default=10000)
    parser.add_argument("--preview-seeds", type=int, default=100)
    parser.add_argument("--jobs", type=int, default=6)
    parser.add_argument("--threads", type=int, default=2)
    args = parser.parse_args()
    if min(args.seeds, args.preview_seeds, args.jobs, args.threads) < 1:
        parser.error("all counts must be positive")
    args.outdir, args.build_dir = args.outdir.resolve(), args.build_dir.resolve()
    history = json.loads(args.history_check.read_text())
    if history.get("passed") is not True:
        raise ValueError("history convergence check has not passed")
    for library in args.build_dir.glob("libgryphon.*"):
        if history["build_sha256"].get(library.name) != digest(library):
            raise ValueError("history check and production must use the same physics build")
    source_paths = [ROOT / "scripts" / p for p in ("run_wide_workflow.py", "run_results.py", "run_halos.py")]
    source_paths += list((ROOT / "plotting/projects/hebreaks").glob("*.py"))
    source_paths += list((ROOT / "plotting/src/gryphon_plots").rglob("*.py"))
    source_paths += list((ROOT / "configs/hebreaks/wide").glob("*.ini"))
    source_hashes = {str(p.relative_to(ROOT)): digest(p) for p in source_paths}
    stages = sorted({min(args.preview_seeds, args.seeds), args.seeds})
    child, stopping = None, False

    def stop(signum, frame):
        nonlocal stopping
        stopping = True
        if child is not None and child.poll() is None:
            child.send_signal(signal.SIGTERM)

    for sig in (signal.SIGINT, signal.SIGTERM):
        signal.signal(sig, stop)

    def status(state, **extra):
        atomic_json(args.outdir / "workflow_status.json", {
            "state": state, "pid": os.getpid(), "updated_utc": datetime.now(timezone.utc).isoformat(),
            "target_per_model": args.seeds, **extra})

    def check_sources():
        for relative, checksum in source_hashes.items():
            if digest(ROOT / relative) != checksum:
                raise RuntimeError(f"workflow source changed: {relative}; refusing mixed analysis")

    with exclusive_run(args.outdir / "workflow_lock") as lock_fd:
        manifest = args.outdir / "workflow_manifest.json"
        identity = {"sources": source_hashes, "history_check": history, "python": sys.version,
                    "build_dir": str(args.build_dir)}
        if manifest.exists() and json.loads(manifest.read_text()) != identity:
            raise RuntimeError("workflow inputs changed; use a fresh output directory")
        atomic_json(manifest, identity)
        try:
            for number in stages:
                if stopping:
                    raise InterruptedError("workflow stopped")
                check_sources()
                status("simulating", current_stage_seeds=number)
                command = [sys.executable, "-u", str(ROOT / "scripts/run_results.py"),
                           "--profile", "wide", "--no-reuse", "--build-dir", str(args.build_dir),
                           "--outdir", str(args.outdir), "--seeds", str(number),
                           "--jobs", str(args.jobs), "--threads", str(args.threads), "--batch-size", "5"]
                child = subprocess.Popen(command, cwd=ROOT, pass_fds=(lock_fd,))
                code = child.wait()
                if code or stopping:
                    raise InterruptedError(f"production stopped (exit {code}); rerun the same workflow to resume")
                check_sources()
                folder = args.outdir / ("analysis" if number == args.seeds else f"analysis_preview_{number}")
                status("analysing", current_stage_seeds=number, analysis_directory=str(folder))
                child = subprocess.Popen([sys.executable, str(ROOT / "plotting/projects/hebreaks/analyse_wide.py"),
                                          "--runs", str(args.outdir), "--seeds", str(number),
                                          "--output", str(folder)], cwd=ROOT, pass_fds=(lock_fd,))
                code = child.wait()
                if stopping:
                    raise InterruptedError("analysis stopped")
                if code:
                    raise RuntimeError(f"analysis failed with exit {code}")
                check_sources()
            status("complete", analysis_directory=str(args.outdir / "analysis"))
            print(f"WORKFLOW COMPLETE: {args.outdir / 'analysis'}", flush=True)
        except InterruptedError as exc:
            status("stopped", error=str(exc))
            print(exc, file=sys.stderr, flush=True)
            return 130
        except Exception as exc:
            status("failed", error=str(exc))
            raise
    return 0


if __name__ == "__main__":
    sys.exit(main())
