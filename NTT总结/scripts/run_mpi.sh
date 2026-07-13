#!/usr/bin/env bash
set -euo pipefail

binary="${1:-./build/ntt_mpi}"
size="${NTT_SIZE:-131072}"
repeats="${NTT_REPEATS:-10}"
threads="${OMP_NUM_THREADS:-2}"

export OMP_NUM_THREADS="$threads"
mpiexec -np 4 "$binary" \
  --size "$size" \
  --threads "$threads" \
  --repeats "$repeats" \
  --target-mod 1337006139375617 \
  --coefficient-limit 1000
