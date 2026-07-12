# main_mpi_stage_v2.cc

保存日期：2026-06-09

版本说明：

- 基于 `main_mpi_crt_modparallel_v1.cc` 继续修改。
- 实现 stage 策略：按 NTT 每一层的 butterfly block 分配给 MPI 进程。
- 每层 butterfly 计算结束后使用 `MPI_Allgatherv` 同步完整数组。
- 点值乘法阶段也按数组区间分配给 MPI 进程，并使用 `MPI_Allgatherv` 同步。
- 为了兼容原有大模数输入，仍保留 CRT 四模数框架；区别是每个 CRT 小模数内部使用 stage MPI 并行，而不是把四个模数分别交给不同进程。
- rank 0 负责最终 CRT 合并、校验、计时输出和写文件。

推荐编译：

```bash
mpic++ -std=c++17 -O3 -fopenmp main.cc -o main
```

推荐运行：

```bash
mpiexec -np 4 -machinefile $PBS_NODEFILE /home/${USER}/main 4 4 5
```
