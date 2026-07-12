#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <pthread.h>
#include <string>
#include <sys/time.h>
#include <unistd.h>
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

int thread_count_for(int jobs) {
    if (jobs <= 1) return 1;

    int threads = 8;
    const char *env = std::getenv("NTT_THREADS");
    if (env != nullptr) {
        int v = std::atoi(env);
        if (v > 0) threads = v;
    } else {
        long cores = sysconf(_SC_NPROCESSORS_ONLN);
        if (cores > 0) threads = (int)cores;
    }

    threads = std::max(1, std::min(threads, 16));
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

void poly_multiply_plain(u64 *a, u64 *b, u64 *ab, int n, u64 mod) {
    std::fill(ab, ab + 2 * n - 1, 0);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            ab[i + j] = (mul_mod(a[i], b[j], mod) + ab[i + j]) % mod;
        }
    }
}

struct NttStageTask {
    u64 *f;
    int len;
    int half;
    int block_l;
    int block_r;
    u64 wn;
    u64 mod;
};

void *ntt_stage_worker(void *arg) {
    NttStageTask *task = (NttStageTask *)arg;
    u64 *f = task->f;
    int len = task->len;
    int half = task->half;
    u64 wn = task->wn;
    u64 mod = task->mod;

    for (int block = task->block_l; block < task->block_r; ++block) {
        int base = block * len;
        u64 w = 1;
        for (int j = 0; j < half; ++j) {
            u64 u = f[base + j];
            u64 v = mul_mod(f[base + j + half], w, mod);
            f[base + j] = add_mod(u, v, mod);
            f[base + j + half] = sub_mod(u, v, mod);
            w = mul_mod(w, wn, mod);
        }
    }
    return nullptr;
}

void run_ntt_stage_pthread(u64 *f, int len, int n, u64 wn, u64 mod) {
    int blocks = n / len;
    int threads = thread_count_for(blocks);

    std::vector<pthread_t> tids(threads);
    std::vector<NttStageTask> tasks(threads);

    for (int t = 0; t < threads; ++t) {
        int l = blocks * t / threads;
        int r = blocks * (t + 1) / threads;
        tasks[t] = {f, len, len >> 1, l, r, wn, mod};
        pthread_create(&tids[t], nullptr, ntt_stage_worker, &tasks[t]);
    }
    for (int t = 0; t < threads; ++t) pthread_join(tids[t], nullptr);
}

struct ScaleTask {
    u64 *f;
    int l;
    int r;
    u64 factor;
    u64 mod;
};

void *scale_worker(void *arg) {
    ScaleTask *task = (ScaleTask *)arg;
    for (int i = task->l; i < task->r; ++i) {
        task->f[i] = mul_mod(task->f[i], task->factor, task->mod);
    }
    return nullptr;
}

void scale_pthread(u64 *f, int n, u64 factor, u64 mod) {
    int threads = thread_count_for(n);
    std::vector<pthread_t> tids(threads);
    std::vector<ScaleTask> tasks(threads);

    for (int t = 0; t < threads; ++t) {
        int l = n * t / threads;
        int r = n * (t + 1) / threads;
        tasks[t] = {f, l, r, factor, mod};
        pthread_create(&tids[t], nullptr, scale_worker, &tasks[t]);
    }
    for (int t = 0; t < threads; ++t) pthread_join(tids[t], nullptr);
}

void ntt_pthread(std::vector<u64> &f, u64 mod, u64 primitive_root, bool inverse) {
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

    for (int len = 2; len <= n; len <<= 1) {
        u64 wn = qpow(primitive_root, (mod - 1) / len, mod);
        if (inverse) wn = qpow(wn, mod - 2, mod);
        run_ntt_stage_pthread(f.data(), len, n, wn, mod);
    }

    if (inverse) {
        u64 inv_n = qpow((u64)n, mod - 2, mod);
        scale_pthread(f.data(), n, inv_n, mod);
    }
}

struct PointwiseTask {
    u64 *a;
    u64 *b;
    int l;
    int r;
    u64 mod;
};

void *pointwise_worker(void *arg) {
    PointwiseTask *task = (PointwiseTask *)arg;
    for (int i = task->l; i < task->r; ++i) {
        task->a[i] = mul_mod(task->a[i], task->b[i], task->mod);
    }
    return nullptr;
}

void pointwise_pthread(std::vector<u64> &a, std::vector<u64> &b, u64 mod) {
    int n = (int)a.size();
    int threads = thread_count_for(n);
    std::vector<pthread_t> tids(threads);
    std::vector<PointwiseTask> tasks(threads);

    for (int t = 0; t < threads; ++t) {
        int l = n * t / threads;
        int r = n * (t + 1) / threads;
        tasks[t] = {a.data(), b.data(), l, r, mod};
        pthread_create(&tids[t], nullptr, pointwise_worker, &tasks[t]);
    }
    for (int t = 0; t < threads; ++t) pthread_join(tids[t], nullptr);
}

void poly_multiply_ntt_pthread(u64 *a, u64 *b, u64 *ab, int n, u64 mod, u64 primitive_root = G) {
    int limit = 1;
    while (limit < 2 * n - 1) limit <<= 1;

    std::vector<u64> A(limit, 0), B(limit, 0);
    for (int i = 0; i < n; ++i) {
        A[i] = a[i] % mod;
        B[i] = b[i] % mod;
    }

    ntt_pthread(A, mod, primitive_root, false);
    ntt_pthread(B, mod, primitive_root, false);
    pointwise_pthread(A, B, mod);
    ntt_pthread(A, mod, primitive_root, true);

    for (int i = 0; i < 2 * n - 1; ++i) ab[i] = A[i] % mod;
}

struct CrtMod {
    u64 mod;
    u64 root;
};

const CrtMod CRT_MODS[4] = {
    {167772161ULL, 3ULL},
    {469762049ULL, 3ULL},
    {754974721ULL, 11ULL},
    {998244353ULL, 3ULL},
};

u64 inv_mod_prime(u64 x, u64 mod) {
    return qpow(x % mod, mod - 2, mod);
}

void prepare_crt_inverse(u64 inv_prod_mod[4]) {
    u128 prod = CRT_MODS[0].mod;
    inv_prod_mod[0] = 1;
    for (int i = 1; i < 4; ++i) {
        u64 mi = CRT_MODS[i].mod;
        inv_prod_mod[i] = inv_mod_prime((u64)(prod % mi), mi);
        prod *= mi;
    }
}

u64 crt_combine_one(const u64 r[4], u64 target_mod, const u64 inv_prod_mod[4]) {
    u128 x = r[0];
    u128 prod = CRT_MODS[0].mod;

    for (int i = 1; i < 4; ++i) {
        u64 mi = CRT_MODS[i].mod;
        u64 cur = (u64)(x % mi);
        u64 need = (r[i] + mi - cur) % mi;
        u64 t = mul_mod(need, inv_prod_mod[i], mi);
        x += prod * t;
        prod *= mi;
    }

    return (u64)(x % target_mod);
}

void poly_multiply_ntt_crt_pthread(u64 *a, u64 *b, u64 *ab, int n, u64 target_mod) {
    int result_len = 2 * n - 1;
    std::vector<std::vector<u64>> residues(4, std::vector<u64>(result_len));
    u64 inv_prod_mod[4];
    prepare_crt_inverse(inv_prod_mod);

    for (int k = 0; k < 4; ++k) {
        poly_multiply_ntt_pthread(a, b, residues[k].data(), n, CRT_MODS[k].mod, CRT_MODS[k].root);
    }

    for (int i = 0; i < result_len; ++i) {
        u64 r[4] = {residues[0][i], residues[1][i], residues[2][i], residues[3][i]};
        ab[i] = crt_combine_one(r, target_mod, inv_prod_mod);
    }
}

void poly_multiply_auto(u64 *a, u64 *b, u64 *ab, int n, u64 mod) {
    if (mod <= 1000000000ULL) {
        poly_multiply_ntt_pthread(a, b, ab, n, mod);
    } else {
        poly_multiply_ntt_crt_pthread(a, b, ab, n, mod);
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
        poly_multiply_auto(a, b, ab, n_, p_);
        auto End = std::chrono::high_resolution_clock::now();

        std::chrono::duration<double, std::milli> elapsed = End - Start;
        fCheck(ab, n_, i);
        std::cout << "latency for n = " << n_ << " p = " << p_
                  << " : " << elapsed.count() << " (ms)" << std::endl;

        fWrite(ab, n_, i);
    }

    return 0;
}
