#!/usr/bin/env bash
set -euo pipefail

binary="${1:-./build/ntt_cpu}"
size="${NTT_SIZE:-131072}"
repeats="${NTT_REPEATS:-10}"
threads="${NTT_THREADS:-4}"

"$binary" --self-test --threads "$threads" >&2
first=1
for backend in serial pthread openmp simd hybrid; do
  if (( first )); then
    "$binary" --backend "$backend" --size "$size" --threads "$threads" --repeats "$repeats"
    first=0
  else
    "$binary" --backend "$backend" --size "$size" --threads "$threads" --repeats "$repeats" | tail -n 1
  fi
done
