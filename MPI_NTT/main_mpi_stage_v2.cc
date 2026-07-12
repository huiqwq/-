#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <mpi.h>
#include <omp.h>
#include <pthread.h>
#include <string>
#include <sys/stat.h>
#include <sys/types.h>
#include <vector>

using u32 = uint32_t;
using u64 = uint64_t;
using u128 = __uint128_t;
using ll = long long;

const u64 G = 3;
const int CRT_CNT = 4;

struct NttMod {
    u64 mod;
    u64 root;
};

const NttMod CRT_MODS[CRT_CNT] = {
    {998244353ULL, 3ULL},
    {1004535809ULL, 3ULL},
    {469762049ULL, 3ULL},
    {1224736769ULL, 3ULL}
};

static inline u64 add_mod(u64 x, u64 y, u64 mod) {
    x += y;
    if (x >= mod || x < y) x -= mod;
    return x;
}

static inline u64 sub_mod(u64 x, u64 y, u64 mod) {
    return x >= y ? x - y : x + mod - y;
}

static inline u64 mul_mod(u64 x, u64 y, u64 mod) {
    return (u64)((u128)x * y % mod);
}

static inline u64 mul_mod_ntt(u64 x, u64 y, u64 mod) {
    return x * y % mod;
}

u64 qpow(u64 a, u64 b, u64 mod) {
    u64 res = 1;
    while (b) {
        if (b & 1) res = mul_mod(res, a, mod);
        a = mul_mod(a, a, mod);
        b >>= 1;
    }
    return res;
}

void fRead(u64 *a, u64 *b, int *n, u64 *p, int input_id) {
    std::string strin = "/nttdata/" + std::to_string(input_id) + ".in";
    std::ifstream fin(strin, std::ios::in);
    fin >> *n >> *p;
    for (int i = 0; i < *n; i++) fin >> a[i];
    for (int i = 0; i < *n; i++) fin >> b[i];
}

void fCheck(u64 *ab, int n, int input_id) {
    std::string strout = "/nttdata/" + std::to_string(input_id) + ".out";
    std::ifstream fin(strout, std::ios::in);
    for (int i = 0; i < n * 2 - 1; i++) {
        u64 x;
        fin >> x;
        if (x != ab[i]) {
            std::cout << "polynomial multiply result is wrong" << std::endl;
            return;
        }
    }
    std::cout << "polynomial multiply result is right" << std::endl;
}

void fWrite(u64 *ab, int n, int input_id) {
    mkdir("files", 0755);
    std::string strout = "files/" + std::to_string(input_id) + ".out";
    std::ofstream fout(strout, std::ios::out);
    for (int i = 0; i < n * 2 - 1; i++) fout << ab[i] << '\n';
}

int get_omp_threads() {
    int threads = 8;
    const char *env = std::getenv("OMP_NUM_THREADS");
    if (env != nullptr) {
        int v = std::atoi(env);
        if (v > 0) threads = v;
    }
    threads = std::max(1, std::min(threads, 8));
    return threads;
}

void init_openmp_runtime() {
    omp_set_dynamic(0);
    omp_set_num_threads(get_omp_threads());
    omp_set_schedule(omp_sched_guided, 64);
}

void bit_reverse(std::vector<u64> &f) {
    int n = (int)f.size();
    for (int i = 1, j = 0; i < n; ++i) {
        int bit = n >> 1;
        while (j & bit) {
            j ^= bit;
            bit >>= 1;
        }
        j ^= bit;

        if (i < j) std::swap(f[i], f[j]);
    }
}

void ntt_openmp(std::vector<u64> &f, u64 mod, u64 root, bool inverse) {
    int n = (int)f.size();
    bit_reverse(f);

    for (int len = 2; len <= n; len <<= 1) {
        u64 wn = qpow(root, (mod - 1) / len, mod);
        if (inverse) wn = qpow(wn, mod - 2, mod);

        int half = len >> 1;
        int blocks = n / len;
        if (blocks >= get_omp_threads()) {
#pragma omp parallel for schedule(static)
            for (int block = 0; block < blocks; ++block) {
                int base = block * len;
                u64 w = 1;
                for (int k = 0; k < half; ++k) {
                    u64 x = f[base + k];
                    u64 y = mul_mod(f[base + k + half], w, mod);
                    f[base + k] = add_mod(x, y, mod);
                    f[base + k + half] = sub_mod(x, y, mod);
                    w = mul_mod(w, wn, mod);
                }
            }
        } else {
            int total = n >> 1;
#pragma omp parallel
            {
                int tid = omp_get_thread_num();
                int nth = omp_get_num_threads();
                int l = total * tid / nth;
                int r = total * (tid + 1) / nth;

                while (l < r) {
                    int block = l / half;
                    int k = l - block * half;
                    int block_end = std::min(r, (block + 1) * half);
                    int base = block * len;
                    u64 w = qpow(wn, (u64)k, mod);

                    for (; l < block_end; ++l, ++k) {
                        u64 x = f[base + k];
                        u64 y = mul_mod(f[base + k + half], w, mod);
                        f[base + k] = add_mod(x, y, mod);
                        f[base + k + half] = sub_mod(x, y, mod);
                        w = mul_mod(w, wn, mod);
                    }
                }
            }
        }
    }

    if (inverse) {
        u64 inv_n = qpow((u64)n, mod - 2, mod);
#pragma omp parallel for schedule(static)
        for (int i = 0; i < n; ++i) {
            f[i] = mul_mod(f[i], inv_n, mod);
        }
    }
}

void ntt_serial(std::vector<u64> &f, u64 mod, u64 root, bool inverse) {
    int n = (int)f.size();
    bit_reverse(f);

    for (int len = 2; len <= n; len <<= 1) {
        u64 wn = qpow(root, (mod - 1) / len, mod);
        if (inverse) wn = qpow(wn, mod - 2, mod);

        int half = len >> 1;
        std::vector<u64> roots(half);
        roots[0] = 1;
        for (int k = 1; k < half; ++k) {
            roots[k] = mul_mod_ntt(roots[k - 1], wn, mod);
        }
        for (int base = 0; base < n; base += len) {
#pragma omp simd
            for (int k = 0; k < half; ++k) {
                u64 x = f[base + k];
                u64 y = mul_mod_ntt(f[base + k + half], roots[k], mod);
                f[base + k] = add_mod(x, y, mod);
                f[base + k + half] = sub_mod(x, y, mod);
            }
        }
    }

    if (inverse) {
        u64 inv_n = qpow((u64)n, mod - 2, mod);
#pragma omp simd
        for (int i = 0; i < n; ++i) {
            f[i] = mul_mod_ntt(f[i], inv_n, mod);
        }
    }
}

void multiply_mod_ntt(const u64 *a, const u64 *b, std::vector<u64> &out, int n, NttMod nm) {
    int limit = 1;
    while (limit < 2 * n - 1) limit <<= 1;

    std::vector<u64> A(limit, 0), B(limit, 0);
    for (int i = 0; i < n; ++i) {
        A[i] = a[i] % nm.mod;
        B[i] = b[i] % nm.mod;
    }

    ntt_serial(A, nm.mod, nm.root, false);
    ntt_serial(B, nm.mod, nm.root, false);
#pragma omp simd
    for (int i = 0; i < limit; ++i) {
        A[i] = mul_mod_ntt(A[i], B[i], nm.mod);
    }
    ntt_serial(A, nm.mod, nm.root, true);

    out.resize(2 * n - 1);
#pragma omp simd
    for (int i = 0; i < 2 * n - 1; ++i) {
        out[i] = A[i];
    }
}

struct ModTask {
    const u64 *a;
    const u64 *b;
    std::vector<u64> *out;
    int n;
    int mod_id;
};

void run_mod_worker(ModTask *task) {
    multiply_mod_ntt(task->a, task->b, *task->out, task->n, CRT_MODS[task->mod_id]);
}

void *mod_worker(void *arg) {
    ModTask *task = (ModTask *)arg;
    run_mod_worker(task);
    return nullptr;
}

void prepare_crt_inverse(u64 inv_prod_mod[CRT_CNT]) {
    u128 prod = CRT_MODS[0].mod;
    inv_prod_mod[0] = 1;
    for (int i = 1; i < CRT_CNT; ++i) {
        inv_prod_mod[i] = qpow((u64)(prod % CRT_MODS[i].mod), CRT_MODS[i].mod - 2, CRT_MODS[i].mod);
        prod *= CRT_MODS[i].mod;
    }
}

u64 crt_combine_one(const u64 r[CRT_CNT], u64 target_mod, const u64 inv_prod_mod[CRT_CNT]) {
    u128 x = r[0];
    u128 prod = CRT_MODS[0].mod;

    for (int i = 1; i < CRT_CNT; ++i) {
        u64 mi = CRT_MODS[i].mod;
        u64 cur = (u64)(x % mi);
        u64 need = (r[i] + mi - cur) % mi;
        u64 t = mul_mod(need, inv_prod_mod[i], mi);
        x += t * prod;
        prod *= mi;
    }
    return (u64)(x % target_mod);
}

struct CombineTask {
    std::vector<std::vector<u64> > *residues;
    u64 *ab;
    int l;
    int r;
    u64 target_mod;
    const u64 *inv_prod_mod;
};

void run_combine_worker(CombineTask *task) {
    std::vector<std::vector<u64> > &residues = *task->residues;
    for (int i = task->l; i < task->r; ++i) {
        u64 r[CRT_CNT];
        for (int j = 0; j < CRT_CNT; ++j) r[j] = residues[j][i];
        task->ab[i] = crt_combine_one(r, task->target_mod, task->inv_prod_mod);
    }
}

void *combine_worker(void *arg) {
    CombineTask *task = (CombineTask *)arg;
    run_combine_worker(task);
    return nullptr;
}

void poly_multiply_crt_openmp_modparallel(u64 *a, u64 *b, u64 *ab, int n, u64 target_mod) {
    int result_len = 2 * n - 1;
    std::vector<std::vector<u64> > residues(CRT_CNT);
    int mod_threads = std::min(get_omp_threads(), CRT_CNT);

#pragma omp parallel for schedule(static) num_threads(mod_threads)
    for (int i = 0; i < CRT_CNT; ++i) {
        multiply_mod_ntt(a, b, residues[i], n, CRT_MODS[i]);
    }

    u64 inv_prod_mod[CRT_CNT];
    prepare_crt_inverse(inv_prod_mod);

#pragma omp parallel for schedule(guided, 256)
    for (int i = 0; i < result_len; ++i) {
        u64 r[CRT_CNT];
        for (int j = 0; j < CRT_CNT; ++j) r[j] = residues[j][i];
        ab[i] = crt_combine_one(r, target_mod, inv_prod_mod);
    }
}

void poly_multiply_crt_openmp(u64 *a, u64 *b, u64 *ab, int n, u64 target_mod) {
    int result_len = 2 * n - 1;
    std::vector<std::vector<u64> > residues(CRT_CNT);

    for (int i = 0; i < CRT_CNT; ++i) {
        multiply_mod_ntt(a, b, residues[i], n, CRT_MODS[i]);
    }

    u64 inv_prod_mod[CRT_CNT];
    prepare_crt_inverse(inv_prod_mod);

#pragma omp parallel for schedule(guided, 256)
    for (int i = 0; i < result_len; ++i) {
        u64 r[CRT_CNT];
        for (int j = 0; j < CRT_CNT; ++j) r[j] = residues[j][i];
        ab[i] = crt_combine_one(r, target_mod, inv_prod_mod);
    }
}

void prepare_counts_by_blocks(int n, int len, int size, std::vector<int> &counts, std::vector<int> &displs) {
    int blocks = n / len;
    counts.assign(size, 0);
    displs.assign(size, 0);

    for (int r = 0; r < size; ++r) {
        int block_l = blocks * r / size;
        int block_r = blocks * (r + 1) / size;
        displs[r] = block_l * len;
        counts[r] = (block_r - block_l) * len;
    }
}

void prepare_counts_by_elements(int n, int size, std::vector<int> &counts, std::vector<int> &displs) {
    counts.assign(size, 0);
    displs.assign(size, 0);

    for (int r = 0; r < size; ++r) {
        int l = n * r / size;
        int rr = n * (r + 1) / size;
        displs[r] = l;
        counts[r] = rr - l;
    }
}

void ntt_mpi_stage(std::vector<u64> &f, u64 mod, u64 root, bool inverse) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int n = (int)f.size();
    bit_reverse(f);

    std::vector<int> counts, displs;
    std::vector<u64> gathered(n);

    for (int len = 2; len <= n; len <<= 1) {
        u64 wn = qpow(root, (mod - 1) / len, mod);
        if (inverse) wn = qpow(wn, mod - 2, mod);

        int half = len >> 1;
        std::vector<u64> roots(half);
        roots[0] = 1;
        for (int k = 1; k < half; ++k) {
            roots[k] = mul_mod_ntt(roots[k - 1], wn, mod);
        }

        int blocks = n / len;
        int block_l = blocks * rank / size;
        int block_r = blocks * (rank + 1) / size;

#pragma omp parallel for schedule(static)
        for (int block = block_l; block < block_r; ++block) {
            int base = block * len;
#pragma omp simd
            for (int k = 0; k < half; ++k) {
                u64 x = f[base + k];
                u64 y = mul_mod_ntt(f[base + k + half], roots[k], mod);
                f[base + k] = add_mod(x, y, mod);
                f[base + k + half] = sub_mod(x, y, mod);
            }
        }

        prepare_counts_by_blocks(n, len, size, counts, displs);
        MPI_Allgatherv(f.data() + displs[rank], counts[rank], MPI_UINT64_T,
                       gathered.data(), counts.data(), displs.data(), MPI_UINT64_T,
                       MPI_COMM_WORLD);
        f.swap(gathered);
    }

    if (inverse) {
        u64 inv_n = qpow((u64)n, mod - 2, mod);
        int l = n * rank / size;
        int r = n * (rank + 1) / size;

#pragma omp parallel for schedule(static)
        for (int i = l; i < r; ++i) {
            f[i] = mul_mod_ntt(f[i], inv_n, mod);
        }

        prepare_counts_by_elements(n, size, counts, displs);
        MPI_Allgatherv(f.data() + displs[rank], counts[rank], MPI_UINT64_T,
                       gathered.data(), counts.data(), displs.data(), MPI_UINT64_T,
                       MPI_COMM_WORLD);
        f.swap(gathered);
    }
}

void pointwise_mpi_stage(std::vector<u64> &A, const std::vector<u64> &B, u64 mod) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int n = (int)A.size();
    int l = n * rank / size;
    int r = n * (rank + 1) / size;

#pragma omp parallel for schedule(static)
    for (int i = l; i < r; ++i) {
        A[i] = mul_mod_ntt(A[i], B[i], mod);
    }

    std::vector<int> counts, displs;
    std::vector<u64> gathered(n);
    prepare_counts_by_elements(n, size, counts, displs);
    MPI_Allgatherv(A.data() + displs[rank], counts[rank], MPI_UINT64_T,
                   gathered.data(), counts.data(), displs.data(), MPI_UINT64_T,
                   MPI_COMM_WORLD);
    A.swap(gathered);
}

void multiply_mod_ntt_mpi_stage(const u64 *a, const u64 *b, std::vector<u64> &out, int n, NttMod nm) {
    int limit = 1;
    while (limit < 2 * n - 1) limit <<= 1;

    std::vector<u64> A(limit, 0), B(limit, 0);
    for (int i = 0; i < n; ++i) {
        A[i] = a[i] % nm.mod;
        B[i] = b[i] % nm.mod;
    }

    ntt_mpi_stage(A, nm.mod, nm.root, false);
    ntt_mpi_stage(B, nm.mod, nm.root, false);
    pointwise_mpi_stage(A, B, nm.mod);
    ntt_mpi_stage(A, nm.mod, nm.root, true);

    out.resize(2 * n - 1);
#pragma omp simd
    for (int i = 0; i < 2 * n - 1; ++i) {
        out[i] = A[i];
    }
}

void poly_multiply_crt_mpi_stage(u64 *a, u64 *b, u64 *ab, int n, u64 target_mod) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int result_len = 2 * n - 1;
    std::vector<std::vector<u64> > residues;
    if (rank == 0) residues.resize(CRT_CNT);

    for (int mod_id = 0; mod_id < CRT_CNT; ++mod_id) {
        std::vector<u64> local;
        multiply_mod_ntt_mpi_stage(a, b, local, n, CRT_MODS[mod_id]);
        if (rank == 0) {
            residues[mod_id].swap(local);
        }
    }

    if (rank == 0) {
        u64 inv_prod_mod[CRT_CNT];
        prepare_crt_inverse(inv_prod_mod);

#pragma omp parallel for schedule(guided, 256)
        for (int i = 0; i < result_len; ++i) {
            u64 r[CRT_CNT];
            for (int j = 0; j < CRT_CNT; ++j) r[j] = residues[j][i];
            ab[i] = crt_combine_one(r, target_mod, inv_prod_mod);
        }
    }
}

void bcast_input_for_mpi(u64 *a, u64 *b, int *n, u64 *p, int rank) {
    MPI_Bcast(n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(p, 1, MPI_UINT64_T, 0, MPI_COMM_WORLD);
    MPI_Bcast(a, *n, MPI_UINT64_T, 0, MPI_COMM_WORLD);
    MPI_Bcast(b, *n, MPI_UINT64_T, 0, MPI_COMM_WORLD);

    if (rank != 0) {
        std::fill(a + *n, a + 300000, 0);
        std::fill(b + *n, b + 300000, 0);
    }
}

void poly_multiply_crt_mpi_modparallel(u64 *a, u64 *b, u64 *ab, int n, u64 target_mod) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int result_len = 2 * n - 1;
    std::vector<std::vector<u64> > residues;
    if (rank == 0) residues.resize(CRT_CNT);

    for (int mod_id = 0; mod_id < CRT_CNT; ++mod_id) {
        int owner = mod_id % size;
        if (rank == owner) {
            std::vector<u64> local;
            multiply_mod_ntt(a, b, local, n, CRT_MODS[mod_id]);

            if (rank == 0) {
                residues[mod_id].swap(local);
            } else {
                MPI_Send(local.data(), result_len, MPI_UINT64_T, 0, 100 + mod_id, MPI_COMM_WORLD);
            }
        }

        if (rank == 0 && owner != 0) {
            residues[mod_id].resize(result_len);
            MPI_Recv(residues[mod_id].data(), result_len, MPI_UINT64_T, owner,
                     100 + mod_id, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    }

    if (rank == 0) {
        u64 inv_prod_mod[CRT_CNT];
        prepare_crt_inverse(inv_prod_mod);

#pragma omp parallel for schedule(guided, 256)
        for (int i = 0; i < result_len; ++i) {
            u64 r[CRT_CNT];
            for (int j = 0; j < CRT_CNT; ++j) r[j] = residues[j][i];
            ab[i] = crt_combine_one(r, target_mod, inv_prod_mod);
        }
    }
}

u64 a[300000], b[300000], ab[300000];

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    init_openmp_runtime();

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int test_begin = 0;
    int test_end = 4;
    int repeat = 1;

    if (argc >= 2) test_begin = test_end = std::atoi(argv[1]);
    if (argc >= 3) {
        test_begin = std::atoi(argv[1]);
        test_end = std::atoi(argv[2]);
    }
    if (argc >= 4) repeat = std::max(1, std::atoi(argv[3]));

    if (rank == 0) {
        std::cout << "MPI stage version, processes = " << size
                  << ", OMP_NUM_THREADS = " << get_omp_threads()
                  << ", repeat = " << repeat << std::endl;
    }

    for (int i = test_begin; i <= test_end; ++i) {
        int n_;
        u64 p_;

        if (rank == 0) {
            fRead(a, b, &n_, &p_, i);
            std::fill(ab, ab + 2 * n_ - 1, 0);
        }

        bcast_input_for_mpi(a, b, &n_, &p_, rank);

        for (int rep = 0; rep < repeat; ++rep) {
            if (rank == 0) std::fill(ab, ab + 2 * n_ - 1, 0);

            MPI_Barrier(MPI_COMM_WORLD);
            double start = MPI_Wtime();
            poly_multiply_crt_mpi_stage(a, b, ab, n_, p_);
            MPI_Barrier(MPI_COMM_WORLD);
            double end = MPI_Wtime();

            if (rank == 0) {
                fCheck(ab, n_, i);
                std::cout << "latency for n = " << n_ << " p = " << p_
                          << " : " << (end - start) * 1000.0 << " (ms)"
                          << " repeat = " << rep << std::endl;
                fWrite(ab, n_, i);
            }
        }
    }

    MPI_Finalize();
    return 0;
}
