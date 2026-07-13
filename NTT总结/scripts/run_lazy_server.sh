#!/usr/bin/env bash
set -euo pipefail

binary="${1:-./build/ntt_lazy_server}"
repeats="${NTT_REPEATS:-10}"
warmups="${NTT_WARMUPS:-2}"
threads="${NTT_THREADS:-8}"

"$binary" --self-test --threads "$threads" >&2
first=1
for n in 1024 16384 131072 262144 1048576; do
  for variant in barrett-serial barrett-lazy-serial barrett-openmp barrett-lazy-openmp; do
    if (( first )); then
      "$binary" --variant "$variant" --size "$n" --threads "$threads" \
        --warmups "$warmups" --repeats "$repeats"
      first=0
    else
      "$binary" --variant "$variant" --size "$n" --threads "$threads" \
        --warmups "$warmups" --repeats "$repeats" | tail -n 1
    fi
  done
done
