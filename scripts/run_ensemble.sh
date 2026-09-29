#!/bin/bash
#
# Run one model over a range of seeds.
#
#   scripts/run_ensemble.sh <executable> <config.ini> <nseeds> [outdir]
#
# e.g. from the build directory:
#
#   ../scripts/run_ensemble.sh ./runYoungPulsars \
#       ../configs/hpwne/youngpulsars_H4.ini 100 ../runs/hpwne

set -euo pipefail

if [ "$#" -lt 3 ]; then
    sed -n '2,10p' "$0" >&2
    exit 1
fi

executable=$1
config=$2
nseeds=$3
outdir=${4:-runs}

for seed in $(seq 0 $((nseeds - 1))); do
    printf '\r%s seed %d/%d' "$(basename "$config")" "$((seed + 1))" "$nseeds"
    "$executable" "$config" --seed "$seed" --outdir "$outdir" >/dev/null
done
printf '\n'
