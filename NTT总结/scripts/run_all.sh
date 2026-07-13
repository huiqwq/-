#!/usr/bin/env bash
set -euo pipefail

root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
mkdir -p "$root/results"

"$root/scripts/run_cpu.sh" "$root/build/ntt_cpu" | tee "$root/results/cpu.csv"
if [[ -x "$root/build/ntt_mpi" ]]; then
  "$root/scripts/run_mpi.sh" "$root/build/ntt_mpi" | tee "$root/results/mpi.csv"
fi
if [[ -x "$root/build/ntt_lazy_server" ]]; then
  "$root/scripts/run_lazy_server.sh" "$root/build/ntt_lazy_server" | tee "$root/results/lazy_server.csv"
fi
if [[ -x "$root/build/hip_ntt_benchmark" ]]; then
  "$root/scripts/run_hip.sh" "$root/build/hip_ntt_benchmark" | tee "$root/results/hip.csv"
fi
