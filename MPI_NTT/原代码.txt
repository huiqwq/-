#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <pthread.h>
#include <string>
#include <vector>
#include <omp.h>
// 可以自行添加需要的头文件
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
    std::string strout = "files/" + std::to_string(input_id) + ".out";
    std::ofstream fout(strout, std::ios::out);
    for (int i = 0; i < n * 2 - 1; i++) fout << ab[i] << '\n';
}
// int get_thread_count(int jobs) {
//     if (jobs <= 1) return 1;

//     int threads = 8;
//     const char *env = std::getenv("NTT_THREADS");
//     if (env != nullptr) {
//         int v = std::atoi(env);
//         if (v > 0) threads = v;
//     }

//     threads = std::max(1, std::min(threads, 8));
//     return std::min(threads, jobs);
// }
int get_omp_threads(){
    int threads = 8;
    const char *env = std::getenv("OMP_NUM_THREADS");
    if (env != nullptr) {
        int v = std::atoi(env);
        if (v > 0) threads = v;
    }
    threads = std::max(1, std::min(threads, 8));
    return threads;
}
void init_openmp_runtime(){
    omp_set_dynamic(1);
    omp_set_num_threads(get_omp_threads());
    omp_set_schedule(omp_sched_guided,64);
}
void bit_reverse(std::vector<u64> &f) {
    int n = (int)f.size();
    for(int i = 1, j = 0; i < n; ++i) {
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
        if(blocks >= get_omp_threads()) {
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

                while(l<r){
                    int block = l / half;
                    int k = l-block*half;
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

void ntt_serial(std::vector<u64> &f, u64 mod, u64 root, bool inverse){
    int n = (int)f.size();
    bit_reverse(f);

    for (int len = 2; len <= n; len <<= 1) {
        u64 wn = qpow(root, (mod - 1) / len, mod);
        if (inverse) wn = qpow(wn, mod - 2, mod);

        int half = len >> 1;
        std::vector<u64> roots(half);//SIMD
        roots[0] = 1;//SIMD
        for (int k = 1; k < half; ++k) {
            roots[k] = mul_mod_ntt(roots[k - 1], wn, mod);
        }//SIMD
        for (int base = 0; base < n; base += len) {
            #pragma omp simd //SIMD
            //u64 w = 1;
            for (int k = 0; k < half; ++k) {
                u64 x = f[base + k];
                //u64 y = mul_mod(f[base + k + half], w, mod);
                u64 y = mul_mod_ntt(f[base + k + half], roots[k], mod);
                f[base + k] = add_mod(x, y, mod);
                f[base + k + half] = sub_mod(x, y, mod);
                //w = mul_mod(w, wn, mod);
            }
        }
    }

    if (inverse) {
        u64 inv_n = qpow((u64)n, mod - 2, mod);
        #pragma omp simd //SIMD
        for (int i = 0; i < n; ++i) {
            //f[i] = mul_mod(f[i], inv_n, mod);
            f[i] = mul_mod_ntt(f[i], inv_n, mod);
        }
    }
}
void multiply_mod_ntt(const u64 *a, const u64 *b, std::vector<u64> &out, int n, NttMod nm){
    int limit = 1;
    while (limit < 2 * n - 1) limit <<= 1;

    std::vector<u64> A(limit, 0), B(limit, 0);
    for (int i = 0; i < n; ++i) {
        A[i] = a[i] % nm.mod;
        B[i] = b[i] % nm.mod;
    }
    // ntt_openmp(A, nm.mod, nm.root, false);
    // ntt_openmp(B, nm.mod, nm.root, false);
    ntt_serial(A, nm.mod, nm.root, false);
    ntt_serial(B, nm.mod, nm.root, false);
    #pragma omp simd//SIMD
    //#pragma omp parallel for schedule(static)
    for (int i = 0; i < limit; ++i) {
        //A[i] = mul_mod(A[i], B[i], nm.mod);
        A[i] = mul_mod_ntt(A[i], B[i], nm.mod);
    }
    //ntt_openmp(A, nm.mod, nm.root, true);
    ntt_serial(A, nm.mod, nm.root, true);
    out.resize(2*n-1);
    #pragma omp simd//SIMD
    for (int i = 0; i < 2 * n - 1; ++i) {
        out[i] = A[i];
    }
}
struct ModTask{
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

u64 crt_combine_one(const u64 r[CRT_CNT], u64 target_mod, const u64 inv_prod_mod[CRT_CNT]){
    u128 x = r[0];
    u128 prod = CRT_MODS[0].mod;

    for(int i=1;i<CRT_CNT;++i){
        u64 mi = CRT_MODS[i].mod;
        u64 cur = (u64)(x%mi);
        u64 need = (r[i]+mi-cur)%mi;
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
// void combine_crt_pthread(std::vector<std::vector<u64> > &residues, u64 *ab, int result_len, u64 target_mod) {
//     u64 inv_prod_mod[CRT_CNT];
//     prepare_crt_inverse(inv_prod_mod);

//     int threads = get_thread_count(result_len);
//     int worker_threads = threads - 1;

//     std::vector<pthread_t> tids(worker_threads);
//     std::vector<CombineTask> tasks(threads);

//     for (int t = 0; t < worker_threads; ++t) {
//         tasks[t].residues = &residues;
//         tasks[t].ab = ab;
//         tasks[t].l = result_len * t / threads;
//         tasks[t].r = result_len * (t + 1) / threads;
//         tasks[t].target_mod = target_mod;
//         tasks[t].inv_prod_mod = inv_prod_mod;
//         pthread_create(&tids[t], nullptr, combine_worker, &tasks[t]);
//     }

//     tasks[worker_threads].residues = &residues;
//     tasks[worker_threads].ab = ab;
//     tasks[worker_threads].l = result_len * worker_threads / threads;
//     tasks[worker_threads].r = result_len;
//     tasks[worker_threads].target_mod = target_mod;
//     tasks[worker_threads].inv_prod_mod = inv_prod_mod;

//     run_combine_worker(&tasks[worker_threads]);

//     for (int t = 0; t < worker_threads; ++t) {
//         pthread_join(tids[t], nullptr);
//     }
// }
// void poly_multiply_crt_pthread(u64 *a, u64 *b, u64 *ab, int n, u64 target_mod) {
//     int result_len = 2 * n - 1;
//     std::vector<std::vector<u64> > residues(CRT_CNT);
//     const int worker_threads = CRT_CNT-1;
//     pthread_t tids[worker_threads];
//     ModTask tasks[CRT_CNT];

//     for (int i = 0; i < worker_threads; ++i) {
//         tasks[i].a = a;
//         tasks[i].b = b;
//         tasks[i].out = &residues[i];
//         tasks[i].n = n;
//         tasks[i].mod_id = i;
//         pthread_create(&tids[i], nullptr, mod_worker, &tasks[i]);
//     }
//     tasks[worker_threads].a = a;
//     tasks[worker_threads].b = b;
//     tasks[worker_threads].out = &residues[worker_threads];
//     tasks[worker_threads].n = n;
//     tasks[worker_threads].mod_id = worker_threads;
//     run_mod_worker(&tasks[worker_threads]);
//     for (int i = 0; i < worker_threads; ++i) {
//         pthread_join(tids[i], nullptr);
//     }

//     combine_crt_pthread(residues, ab, result_len, target_mod);
// }
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

    #pragma omp parallel for schedule(guided,256)
    for (int i = 0; i < result_len; ++i) {
        u64 r[CRT_CNT];
        for (int j = 0; j < CRT_CNT; ++j) r[j] = residues[j][i];
        ab[i] = crt_combine_one(r, target_mod, inv_prod_mod);
    }
}
// void poly_multiply_crt_serial(u64 *a, u64 *b, u64 *ab, int n, u64 target_mod) {
//     int result_len = 2 * n - 1;
//     std::vector<std::vector<u64> > residues(CRT_CNT);

//     for (int i = 0; i < CRT_CNT; ++i) {
//         multiply_mod_ntt(a, b, residues[i], n, CRT_MODS[i]);
//     }

//     u64 inv_prod_mod[CRT_CNT];
//     prepare_crt_inverse(inv_prod_mod);

//     for (int i = 0; i < result_len; ++i) {
//         u64 r[CRT_CNT];
//         for (int j = 0; j < CRT_CNT; ++j) r[j] = residues[j][i];
//         ab[i] = crt_combine_one(r, target_mod, inv_prod_mod);
//     }
// }
// struct NttContext {
//     u64 *f;
//     int n;
//     int threads;
//     u64 mod;
//     bool inverse;
//     pthread_barrier_t barrier;
// };

// struct NttTask {
//     NttContext *ctx;
//     int tid;
// };

// void ntt_stage_by_blocks(u64 *f, int n, int len, u64 wn, u64 mod, int tid, int threads) {
//     int half = len >> 1;
//     int blocks = n / len;

//     for (int block = tid; block < blocks; block += threads) {
//         int base = block * len;
//         u64 w = 1;
//         for (int k = 0; k < half; ++k) {
//             u64 x = f[base + k];
//             u64 y = mul_mod(f[base + k + half], w, mod);
//             f[base + k] = add_mod(x, y, mod);
//             f[base + k + half] = sub_mod(x, y, mod);
//             w = mul_mod(w, wn, mod);
//         }
//     }
// }

// void ntt_stage_by_butterflies(u64 *f, int n, int len, u64 wn, u64 mod, int tid, int threads) {
//     int half = len >> 1;
//     int total = n >> 1;
//     int l = total * tid / threads;
//     int r = total * (tid + 1) / threads;

//     while (l < r) {
//         int block = l / half;
//         int k = l % half;
//         int block_end = std::min(r, (block + 1) * half);
//         int base = block * len;
//         u64 w = qpow(wn, k, mod);

//         for (; l < block_end; ++l, ++k) {
//             u64 x = f[base + k];
//             u64 y = mul_mod(f[base + k + half], w, mod);
//             f[base + k] = add_mod(x, y, mod);
//             f[base + k + half] = sub_mod(x, y, mod);
//             w = mul_mod(w, wn, mod);
//         }
//     }
// }

// void run_ntt_worker(NttContext *ctx, int tid) {
//     u64 *f = ctx->f;
//     int n = ctx->n;
//     int threads = ctx->threads;
//     u64 mod = ctx->mod;

//     for (int len = 2; len <= n; len <<= 1) {
//         u64 wn = qpow(G, (mod - 1) / len, mod);
//         if (ctx->inverse) wn = qpow(wn, mod - 2, mod);

//         int blocks = n / len;
//         if (blocks >= threads) {
//             ntt_stage_by_blocks(f, n, len, wn, mod, tid, threads);
//         } else {
//             ntt_stage_by_butterflies(f, n, len, wn, mod, tid, threads);
//         }
//         pthread_barrier_wait(&ctx->barrier);
//     }

//     if (ctx->inverse) {
//         u64 inv_n = qpow((u64)n, mod - 2, mod);
//         for (int i = tid; i < n; i += threads) {
//             f[i] = mul_mod(f[i], inv_n, mod);
//         }
//         pthread_barrier_wait(&ctx->barrier);
//     }
// }

// void *ntt_worker(void *arg) {
//     NttTask *task = (NttTask *)arg;
//     run_ntt_worker(task->ctx, task->tid);
//     return nullptr;
// }

// void ntt_pthread(std::vector<u64> &f, u64 mod, bool inverse) {
//     int n = (int)f.size();
//     bit_reverse(f);

//     int threads = get_thread_count(n >> 1);
//     int worker_threads = threads - 1;
//     NttContext ctx;
//     ctx.f = f.data();
//     ctx.n = n;
//     ctx.threads = threads;
//     ctx.mod = mod;
//     ctx.inverse = inverse;
//     pthread_barrier_init(&ctx.barrier, nullptr, threads);

//     std::vector<pthread_t> tids(worker_threads);
//     std::vector<NttTask> tasks(worker_threads);

//     for (int t = 0; t < worker_threads; ++t) {
//         tasks[t] = {&ctx, t};
//         pthread_create(&tids[t], nullptr, ntt_worker, &tasks[t]);
//     }

//     run_ntt_worker(&ctx, worker_threads);

//     for (int t = 0; t < worker_threads; ++t) {
//         pthread_join(tids[t], nullptr);
//     }

//     pthread_barrier_destroy(&ctx.barrier);
// }

// struct PointwiseTask {
//     u64 *a;
//     u64 *b;
//     int n;
//     int tid;
//     int threads;
//     u64 mod;
// };

// void run_pointwise_worker(PointwiseTask *task) {
//     for (int i = task->tid; i < task->n; i += task->threads) {
//         task->a[i] = mul_mod(task->a[i], task->b[i], task->mod);
//     }
// }

// void *pointwise_worker(void *arg) {
//     PointwiseTask *task = (PointwiseTask *)arg;
//     run_pointwise_worker(task);
//     return nullptr;
// }

// void pointwise_pthread(std::vector<u64> &a, std::vector<u64> &b, u64 mod) {
//     int n = (int)a.size();
//     int threads = get_thread_count(n);
//     int worker_threads = threads - 1;

//     std::vector<pthread_t> tids(worker_threads);
//     std::vector<PointwiseTask> tasks(threads);

//     for (int t = 0; t < worker_threads; ++t) {
//         tasks[t] = {a.data(), b.data(), n, t, threads, mod};
//         pthread_create(&tids[t], nullptr, pointwise_worker, &tasks[t]);
//     }

//     tasks[worker_threads] = {a.data(), b.data(), n, worker_threads, threads, mod};
//     run_pointwise_worker(&tasks[worker_threads]);

//     for (int t = 0; t < worker_threads; ++t) {
//         pthread_join(tids[t], nullptr);
//     }
// }

// void poly_multiply_ntt_pthread(u64 *a, u64 *b, u64 *ab, int n, u64 mod) {
//     int limit = 1;
//     while (limit < 2 * n - 1) limit <<= 1;

//     std::vector<u64> A(limit, 0), B(limit, 0);
//     for (int i = 0; i < n; ++i) {
//         A[i] = a[i] % mod;
//         B[i] = b[i] % mod;
//     }

//     ntt_pthread(A, mod, false);
//     ntt_pthread(B, mod, false);
//     pointwise_pthread(A, B, mod);
//     ntt_pthread(A, mod, true);

//     for (int i = 0; i < 2 * n - 1; ++i) {
//         ab[i] = A[i] % mod;
//     }
// }


u64 a[300000], b[300000], ab[300000];
int main(int argc, char *argv[]) {
    int test_begin = 0;
    int test_end = 4;
    // 保证输入的所有模数的原根均为 3, 且模数都能表示为 a \times 4 ^ k + 1 的形式
    // 输入模数分别为 7340033 104857601 469762049 263882790666241
    // 第四个模数超过了整型表示范围, 如果实现此模数意义下的多项式乘法需要修改框架
    // 对第四个模数的输入数据不做必要要求, 如果要自行探索大模数 NTT, 请在完成前三个模数的基础代码及优化后实现大模数 NTT
    // 输入文件共五个, 第一个输入文件 n = 4, 其余四个文件分别对应四个模数, n = 131072
    // 在实现快速数论变化前, 后四个测试样例运行时间较久, 推荐调试正确性时只使用输入文件 1
    

    for (int i = test_begin; i <= test_end; ++i) {
        int n_;
        u64 p_;

        fRead(a, b, &n_, &p_, i);
        std::fill(ab, ab + 2 * n_ - 1, 0);

        auto Start = std::chrono::high_resolution_clock::now();
        poly_multiply_crt_openmp_modparallel(a, b, ab, n_, p_);
        auto End = std::chrono::high_resolution_clock::now();

        std::chrono::duration<double, std::milli> elapsed = End - Start;
        fCheck(ab, n_, i);
        std::cout << "latency for n = " << n_ << " p = " << p_
                  << " : " << elapsed.count() << " (ms)" << std::endl;

        fWrite(ab, n_, i);
    }
    return 0;
}
