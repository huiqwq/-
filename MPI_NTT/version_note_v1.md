# main_mpi_crt_modparallel_v1.cc

保存日期：2026-06-09

版本说明：

- 基于上次实验的 OpenMP + SIMD + CRT NTT 代码框架修改。
- 保留原有 `fRead/fCheck/fWrite`、`ntt_serial`、`multiply_mod_ntt`、CRT 合并和主测试循环。
- 新增 MPI 初始化、输入广播、CRT 四模数任务分配、结果收集。
- 并行策略：四个 CRT 模数按 `mod_id % MPI进程数` 分配给不同进程。
- rank 0 负责读取输入、CRT 合并、校验、计时输出、写文件。
- 推荐编译：

```bash
mpic++ -std=c++17 -O3 -fopenmp main.cc -o main
```

- 推荐运行：

```bash
mpiexec -np 4 -machinefile $PBS_NODEFILE /home/${USER}/main 0 4 1
```
