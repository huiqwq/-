#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <omp.h>
#include <string>
#include <vector>

using u64 = uint64_t;
using u128 = __uint128_t;

const int MAXN = 300000;

u64 a[MAXN], b[MAXN], ab[MAXN];

int get_omp_threads() {
    int threads = 8;
    const char *env = std::getenv("OMP_NUM_THREADS");
    if (env != nullptr) {
        int v = std::atoi(env);
        if (v > 0) threads = v;
    }
    return std::max(1, std::min(threads, 8));
}

void init_openmp_runtime() {
    omp_set_dynamic(0);
    omp_set_num_threads(get_omp_threads());
    omp_set_schedule(omp_sched_guided, 256);
}

struct NttMod {
    u64 mod;
    u64 root;
};

const int CRT_CNT = 4;
const NttMod CRT_MODS[CRT_CNT] = {
    {998244353ULL, 3ULL},
    {1004535809ULL, 3ULL},
    {469762049ULL, 3ULL},
    {1224736769ULL, 3ULL},
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

void ntt_serial(std::vector<u64> &f, u64 mod, u64 root, bool inverse) {
    int n = (int)f.size();
    bit_reverse(f);

    for (int len = 2; len <= n; len <<= 1) {
        u64 wn = qpow(root, (mod - 1) / len, mod);
        if (inverse) wn = qpow(wn, mod - 2, mod);

        int half = len >> 1;
        for (int base = 0; base < n; base += len) {
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

    if (inverse) {
        u64 inv_n = qpow((u64)n, mod - 2, mod);
        for (int i = 0; i < n; ++i) {
            f[i] = mul_mod(f[i], inv_n, mod);
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
    for (int i = 0; i < limit; ++i) {
        A[i] = mul_mod(A[i], B[i], nm.mod);
    }
    ntt_serial(A, nm.mod, nm.root, true);

    int result_len = 2 * n - 1;
    out.resize(result_len);
    for (int i = 0; i < result_len; ++i) {
        out[i] = A[i];
    }
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
        x += prod * t;
        prod *= mi;
    }

    return (u64)(x % target_mod);
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

int main(int argc, char *argv[]) {
    init_openmp_runtime();

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
