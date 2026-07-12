# main_mpi_stage_barrett_difdit_v4.cc

保存日期：2026-06-09

版本说明：

- 基于 `main_mpi_stage_barrett_v3.cc` 继续修改。
- 保留 MPI stage 策略：NTT 每层 butterfly block 按 MPI 进程划分，每层后 `MPI_Allgatherv` 同步。
- 保留 Barrett 规约优化，用于 CRT 小模数 NTT 内部的高频模乘。
- 新增 DIF/DIT 优化：
  - 正变换使用 DIF，从 `len = n` 递减到 `2`。
  - DIF 正变换不做初始 `bit_reverse`，输出为 bit-reversed 顺序。
  - 两个输入多项式都以相同 bit-reversed 顺序输出，因此可以直接逐点相乘。
  - 逆变换使用 DIT，从 `len = 2` 递增到 `n`，直接接收 bit-reversed 输入并输出正常顺序结果。
  - 因此主 NTT 流程避免了显式位逆序置换。
- CRT 大模数合并仍使用原来的 `__uint128_t` 逻辑。

推荐编译：

```bash
mpic++ -std=c++17 -O3 -fopenmp main.cc -o main
```

推荐运行：

```bash
mpiexec -np 4 -machinefile $PBS_NODEFILE /home/${USER}/main 4 4 5
```
