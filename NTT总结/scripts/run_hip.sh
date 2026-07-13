#!/usr/bin/env bash
set -euo pipefail

root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
binary="${1:-$root/build/hip_ntt_benchmark}"
repeats="${NTT_REPEATS:-10}"
warmups="${NTT_WARMUPS:-2}"

mkdir -p "$root/build"
if [[ ! -x "$binary" ]]; then
  hipcc -O3 -std=c++17 -I"$root/include" "$root/src/hip_ntt_benchmark.cpp" -o "$binary"
fi

"$binary" --self-test --block-size 128 >&2
first=1
for n in 1024 16384 131072 262144 1048576; do
  for variant in runtime barrett montgomery tiled lazy; do
    if (( first )); then
      "$binary" --variant "$variant" --size "$n" --block-size 128 \
        --warmups "$warmups" --repeats "$repeats"
      first=0
    else
      "$binary" --variant "$variant" --size "$n" --block-size 128 \
        --warmups "$warmups" --repeats "$repeats" | tail -n 1
    fi
  done
done
