#!/usr/bin/env bash
# One independent job per structure (runs in parallel on as many nodes as the queue allows).
#   bash submit_all.sh              -> runs/*       (dry ribbons: DZVP optimisation + charges)
#   bash submit_all.sh runs_wet     -> runs_wet/*   (ribbon + water: single points)
# Extra qsub options can be appended, e.g.  bash submit_all.sh runs "-l walltime=48:00:00".
cd "$(dirname "${BASH_SOURCE[0]}")"
dir="${1:-runs}"; shift || true
for d in "$dir"/*/; do
  n=$(basename "$d")
  [[ -f "$d/charges.dat" ]] && { echo "skip $n (done)"; continue; }
  qsub -N "c2_${n:0:12}" -v NAME="$n",RUNS="$dir" "$@" run_cp2k_serial.pbs
done
