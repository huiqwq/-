#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <pthread.h>
#include <string>
#include <vector>

using u64 = uint64_t;
using u128 = __uint128_t;

const int MAXN = 300000;
const u64 G = 3;

u64 a[MAXN], b[MAXN], ab[MAXN];

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

u64 qpow(u64 a, u64 b, u64 mod) {
    u64 res = 1;
    while (b) {
        if (b & 1) res = mul_mod(res, a, mod);
        a = mul_mod(a, a, mod);
        b >>= 1;
    }
    return res;
}

int get_thread_count(int jobs) {
    if (jobs <= 1) return 1;

    int threads = 8;
    const char *env = std::getenv("NTT_THREADS");
    if (env != nullptr) {
        int v = std::atoi(env);
        if (v > 0) threads = v;
    }

    threads = std::max(1, std::min(threads, 8));
    return std::min(threads, jobs);
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

struct NttContext {
    u64 *f;
    int n;
    int threads;
    u64 mod;
    bool inverse;
    pthread_barrier_t barrier;
};

struct NttTask {
    NttContext *ctx;
    int tid;
};

void ntt_stage_by_blocks(u64 *f, int n, int len, u64 wn, u64 mod, int tid, int threads) {
    int half = len >> 1;
    int blocks = n / len;

    for (int block = tid; block < blocks; block += threads) {
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
}

void ntt_stage_by_butterflies(u64 *f, int n, int len, u64 wn, u64 mod, int tid, int threads) {
    int half = len >> 1;
    int total = n >> 1;
    int l = total * tid / threads;
    int r = total * (tid + 1) / threads;

    while (l < r) {
        int block = l / half;
        int k = l % half;
        int block_end = std::min(r, (block + 1) * half);
        int base = block * len;
        u64 w = qpow(wn, k, mod);

        for (; l < block_end; ++l, ++k) {
            u64 x = f[base + k];
            u64 y = mul_mod(f[base + k + half], w, mod);
            f[base + k] = add_mod(x, y, mod);
            f[base + k + half] = sub_mod(x, y, mod);
            w = mul_mod(w, wn, mod);
        }
    }
}

void run_ntt_worker(NttContext *ctx, int tid) {
    u64 *f = ctx->f;
    int n = ctx->n;
    int threads = ctx->threads;
    u64 mod = ctx->mod;

    for (int len = 2; len <= n; len <<= 1) {
        u64 wn = qpow(G, (mod - 1) / len, mod);
        if (ctx->inverse) wn = qpow(wn, mod - 2, mod);

        int blocks = n / len;
        if (blocks >= threads) {
            ntt_stage_by_blocks(f, n, len, wn, mod, tid, threads);
        } else {
            ntt_stage_by_butterflies(f, n, len, wn, mod, tid, threads);
        }
        pthread_barrier_wait(&ctx->barrier);
    }

    if (ctx->inverse) {
        u64 inv_n = qpow((u64)n, mod - 2, mod);
        for (int i = tid; i < n; i += threads) {
            f[i] = mul_mod(f[i], inv_n, mod);
        }
        pthread_barrier_wait(&ctx->barrier);
    }
}

void *ntt_worker(void *arg) {
    NttTask *task = (NttTask *)arg;
    run_ntt_worker(task->ctx, task->tid);
    return nullptr;
}

void ntt_pthread(std::vector<u64> &f, u64 mod, bool inverse) {
    int n = (int)f.size();
    bit_reverse(f);

    int threads = get_thread_count(n >> 1);
    int worker_threads = threads - 1;
    NttContext ctx;
    ctx.f = f.data();
    ctx.n = n;
    ctx.threads = threads;
    ctx.mod = mod;
    ctx.inverse = inverse;
    pthread_barrier_init(&ctx.barrier, nullptr, threads);

    std::vector<pthread_t> tids(worker_threads);
    std::vector<NttTask> tasks(worker_threads);

    for (int t = 0; t < worker_threads; ++t) {
        tasks[t] = {&ctx, t};
        pthread_create(&tids[t], nullptr, ntt_worker, &tasks[t]);
    }

    run_ntt_worker(&ctx, worker_threads);

    for (int t = 0; t < worker_threads; ++t) {
        pthread_join(tids[t], nullptr);
    }

    pthread_barrier_destroy(&ctx.barrier);
}

struct PointwiseTask {
    u64 *a;
    u64 *b;
    int n;
    int tid;
    int threads;
    u64 mod;
};

void run_pointwise_worker(PointwiseTask *task) {
    for (int i = task->tid; i < task->n; i += task->threads) {
        task->a[i] = mul_mod(task->a[i], task->b[i], task->mod);
    }
}

void *pointwise_worker(void *arg) {
    PointwiseTask *task = (PointwiseTask *)arg;
    run_pointwise_worker(task);
    return nullptr;
}

void pointwise_pthread(std::vector<u64> &a, std::vector<u64> &b, u64 mod) {
    int n = (int)a.size();
    int threads = get_thread_count(n);
    int worker_threads = threads - 1;

    std::vector<pthread_t> tids(worker_threads);
    std::vector<PointwiseTask> tasks(threads);

    for (int t = 0; t < worker_threads; ++t) {
        tasks[t] = {a.data(), b.data(), n, t, threads, mod};
        pthread_create(&tids[t], nullptr, pointwise_worker, &tasks[t]);
    }

    tasks[worker_threads] = {a.data(), b.data(), n, worker_threads, threads, mod};
    run_pointwise_worker(&tasks[worker_threads]);

    for (int t = 0; t < worker_threads; ++t) {
        pthread_join(tids[t], nullptr);
    }
}

void poly_multiply_ntt_pthread(u64 *a, u64 *b, u64 *ab, int n, u64 mod) {
    int limit = 1;
    while (limit < 2 * n - 1) limit <<= 1;

    std::vector<u64> A(limit, 0), B(limit, 0);
    for (int i = 0; i < n; ++i) {
        A[i] = a[i] % mod;
        B[i] = b[i] % mod;
    }

    ntt_pthread(A, mod, false);
    ntt_pthread(B, mod, false);
    pointwise_pthread(A, B, mod);
    ntt_pthread(A, mod, true);

    for (int i = 0; i < 2 * n - 1; ++i) {
        ab[i] = A[i] % mod;
    }
}

int main(int argc, char *argv[]) {
    int test_begin = 0;
    int test_end = 4;

    if (argc >= 2) test_begin = test_end = std::atoi(argv[1]);
    if (argc >= 3) {
        test_begin = std::atoi(argv[1]);
        test_end = std::atoi(argv[2]);
    }

    for (int i = test_begin; i <= test_end; ++i) {
        int n_;
        u64 p_;

        fRead(a, b, &n_, &p_, i);
        std::fill(ab, ab + 2 * n_ - 1, 0);

        auto Start = std::chrono::high_resolution_clock::now();
        poly_multiply_ntt_pthread(a, b, ab, n_, p_);
        auto End = std::chrono::high_resolution_clock::now();

        std::chrono::duration<double, std::milli> elapsed = End - Start;
        fCheck(ab, n_, i);
        std::cout << "latency for n = " << n_ << " p = " << p_
                  << " : " << elapsed.count() << " (ms)" << std::endl;

        fWrite(ab, n_, i);
    }

    return 0;
}
