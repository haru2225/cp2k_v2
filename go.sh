#!/usr/bin/env bash
# Submit with more nodes per structure/system:   bash go.sh [nodes=4] [walltime=24:00:00] [extra qsub options]
# Each sub-job then runs CP2K on nodes x 16 MPI ranks (sc16 = 16 cores per node). More nodes = faster per structure,
# but a bigger request waits longer in the queue; 2-4 nodes is a reasonable range for 90-150 atoms.
cd "$(dirname "${BASH_SOURCE[0]}")"
nodes="${1:-4}"; wt="${2:-24:00:00}"; [[ $# -ge 2 ]] && shift 2 || shift $#
qsub -l select=${nodes}:ncpus=16:mpiprocs=16:ompthreads=1 -l walltime="$wt" "$@" run_cp2k.pbs
