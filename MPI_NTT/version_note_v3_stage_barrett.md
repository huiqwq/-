# main_mpi_stage_barrett_v3.cc

保存日期：2026-06-09

版本说明：

- 基于 `main_mpi_stage_v2.cc` 继续修改。
- 保留 MPI stage 策略：NTT 每层 butterfly block 按 MPI 进程划分，每层后 `MPI_Allgatherv` 同步。
- 新增 Barrett 规约优化，用于 CRT 小模数 NTT 内部的高频模乘。
- Barrett 参数：
  - `mu = floor(2^64 / mod)`
  - `q = high64(x * mu)`
  - `r = x - q * mod`
  - 使用两次条件减法修正到 `[0, mod)`。
- 替换位置：
  - 旋转因子 `roots[k]` 递推。
  - butterfly 中的 `f[i] * roots[k] mod mod`。
  - 点值乘法 `A[i] * B[i] mod mod`。
  - 逆 NTT 中乘 `inv_n` 的缩放。
- CRT 大模数合并仍使用原来的 `__uint128_t` 逻辑，未使用 Barrett。

推荐编译：

```bash
mpic++ -std=c++17 -O3 -fopenmp main.cc -o main
```

推荐运行：

```bash
mpiexec -np 4 -machinefile $PBS_NODEFILE /home/${USER}/main 4 4 5
```
