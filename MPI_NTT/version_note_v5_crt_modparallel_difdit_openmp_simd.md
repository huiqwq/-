# main_mpi_crt_modparallel_difdit_openmp_simd_v5.cc

保存日期：2026-06-09

版本说明：

- 基于 `main_mpi_crt_modparallel_v1.cc` 继续修改。
- 保留 v1 的 CRT 模数级 MPI 并行策略：
  - 四个 CRT 小模数任务按 `mod_id % MPI进程数` 分配。
  - 非 0 进程计算 residue 后发送给 rank 0。
  - rank 0 负责 CRT 合并、校验、计时输出和写文件。
- 新增进程内 OpenMP + SIMD 优化：
  - 每个 MPI 进程处理自己的小模数 NTT 时，内部使用 `#pragma omp parallel for`。
  - butterfly 内层循环、点值乘法、结果拷贝使用 `#pragma omp simd` 或 `parallel for simd`。
- 新增 DIF/DIT 优化：
  - 正变换使用 DIF，不做初始 `bit_reverse`。
  - DIF 正变换输出 bit-reversed 顺序。
  - 两个输入数组的顺序一致，可以直接逐点相乘。
  - 逆变换使用 DIT，接收 bit-reversed 输入并输出正常顺序结果。
- 与 v4 的区别：
  - v4 是 MPI stage 并行，每层使用 `MPI_Allgatherv`。
  - v5 是 v1 的低通信 CRT 模数级 MPI 并行，仅在每个进程内部使用 OpenMP/SIMD 和 DIF/DIT。

推荐编译：

```bash
mpic++ -std=c++17 -O3 -fopenmp main.cc -o main
```

推荐运行：

```bash
mpiexec -np 4 -machinefile $PBS_NODEFILE /home/${USER}/main 4 4 5
```
