#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <pthread.h>
#include <string>
#include <vector>

using u64 = uint64_t;
using u128 = __uint128_t;

const int MAXN = 300000;

u64 a[MAXN], b[MAXN], ab[MAXN];

static inline u64 mul_mod(u64 x, u64 y, u64 mod) {
    return (u64)((u128)x * y % mod);
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

struct PlainTask {
    const u64 *a;
    const u64 *b;
    u64 *ab;
    int n;
    int l;
    int r;
    u64 mod;
};

void *plain_worker(void *arg) {
    PlainTask *task = (PlainTask *)arg;
    const u64 *a = task->a;
    const u64 *b = task->b;
    u64 *ab = task->ab;
    int n = task->n;
    u64 mod = task->mod;

    for (int k = task->l; k < task->r; ++k) {
        int begin = std::max(0, k - n + 1);
        int end = std::min(n - 1, k);

        u128 sum = 0;
        for (int i = begin; i <= end; ++i) {
            sum += (u128)a[i] * b[k - i];
            if ((i & 255) == 255) sum %= mod;
        }
        ab[k] = (u64)(sum % mod);
    }

    return nullptr;
}

void poly_multiply_plain_pthread(u64 *a, u64 *b, u64 *ab, int n, u64 mod) {
    int result_len = 2 * n - 1;
    int threads = get_thread_count(result_len);
    int worker_threads = threads - 1;

    std::vector<pthread_t> tids(worker_threads);
    std::vector<PlainTask> tasks(threads);

    for (int t = 0; t < worker_threads; ++t) {
        int l = result_len * t / threads;
        int r = result_len * (t + 1) / threads;
        tasks[t] = {a, b, ab, n, l, r, mod};
        pthread_create(&tids[t], nullptr, plain_worker, &tasks[t]);
    }

    int main_l = result_len * worker_threads / threads;
    int main_r = result_len;
    tasks[worker_threads] = {a, b, ab, n, main_l, main_r, mod};
    plain_worker(&tasks[worker_threads]);

    for (int t = 0; t < worker_threads; ++t) {
        pthread_join(tids[t], nullptr);
    }
}

void poly_multiply_plain_serial(u64 *a, u64 *b, u64 *ab, int n, u64 mod) {
    std::fill(ab, ab + 2 * n - 1, 0);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            ab[i + j] = (ab[i + j] + mul_mod(a[i], b[j], mod)) % mod;
        }
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
        poly_multiply_plain_pthread(a, b, ab, n_, p_);
        auto End = std::chrono::high_resolution_clock::now();

        std::chrono::duration<double, std::milli> elapsed = End - Start;
        fCheck(ab, n_, i);
        std::cout << "latency for n = " << n_ << " p = " << p_
                  << " : " << elapsed.count() << " (ms)" << std::endl;

        fWrite(ab, n_, i);
    }

    return 0;
}
