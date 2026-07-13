#include "ntt_common.hpp"

#include <pthread.h>

#include <cstdlib>
#include <cstring>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <memory>
#include <numeric>
#include <string_view>

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef __AVX2__
#include <immintrin.h>
#endif

namespace {

using ntt::Plan;
using ntt::u32;

enum class Backend { Serial, Pthread, OpenMP, Simd, Hybrid };

struct Options {
    std::size_t transform_size = 1u << 17u;
    int threads = 4;
    int repeats = 5;
    Backend backend = Backend::Serial;
    bool self_test = false;
};

const char* backend_name(Backend backend) {
    switch (backend) {
        case Backend::Serial: return "serial";
        case Backend::Pthread: return "pthread";
        case Backend::OpenMP: return "openmp";
        case Backend::Simd: return "simd";
        case Backend::Hybrid: return "openmp+simd";
    }
    return "unknown";
}

Backend parse_backend(std::string_view value) {
    if (value == "serial") return Backend::Serial;
    if (value == "pthread") return Backend::Pthread;
    if (value == "openmp") return Backend::OpenMP;
    if (value == "simd") return Backend::Simd;
    if (value == "hybrid") return Backend::Hybrid;
    throw std::invalid_argument("backend must be serial, pthread, openmp, simd, or hybrid");
}

Options parse_options(int argc, char** argv) {
    Options options;
    for (int i = 1; i < argc; ++i) {
        const std::string_view arg(argv[i]);
        auto require_value = [&](const char* name) -> const char* {
            if (++i == argc) throw std::invalid_argument(std::string("missing value for ") + name);
            return argv[i];
        };
        if (arg == "--size") options.transform_size = std::stoull(require_value("--size"));
        else if (arg == "--threads") options.threads = std::stoi(require_value("--threads"));
        else if (arg == "--repeats") options.repeats = std::stoi(require_value("--repeats"));
        else if (arg == "--backend") options.backend = parse_backend(require_value("--backend"));
        else if (arg == "--self-test") options.self_test = true;
        else if (arg == "--help") {
            std::cout << "Usage: ntt_cpu [--backend serial|pthread|openmp|simd|hybrid] "
                         "[--size N] [--threads T] [--repeats R] [--self-test]\n";
            std::exit(0);
        } else {
            throw std::invalid_argument("unknown argument: " + std::string(arg));
        }
    }
    if (!ntt::is_power_of_two(options.transform_size) || options.transform_size < 2) {
        throw std::invalid_argument("--size must be a power of two >= 2");
    }
    if (options.threads < 1 || options.repeats < 1) {
        throw std::invalid_argument("threads and repeats must be positive");
    }
    return options;
}

struct PthreadContext {
    std::vector<u32>* data{};
    const Plan* plan{};
    bool inverse{};
    int thread_count{};
    pthread_barrier_t barrier{};
};

struct WorkerArgs {
    PthreadContext* context{};
    int id{};
};

void* pthread_worker(void* raw) {
    auto* args = static_cast<WorkerArgs*>(raw);
    PthreadContext& ctx = *args->context;
    auto& a = *ctx.data;
    const auto& plan = *ctx.plan;
    const u32 mod = plan.modulus.value;
    const auto& stages = ctx.inverse ? plan.inverse_twiddles : plan.forward_twiddles;

    for (std::size_t stage = 0, len = 2; stage < stages.size(); ++stage, len <<= 1u) {
        const std::size_t butterflies = plan.n / 2;
        const std::size_t begin = butterflies * args->id / ctx.thread_count;
        const std::size_t end = butterflies * (args->id + 1) / ctx.thread_count;
        const auto& twiddle = stages[stage];
        for (std::size_t butterfly = begin; butterfly < end; ++butterfly) {
            const std::size_t j = butterfly % (len / 2);
            const std::size_t base = (butterfly / (len / 2)) * len;
            const u32 u = a[base + j];
            const u32 v = ntt::mul_mod(a[base + j + len / 2], twiddle[j], mod);
            a[base + j] = ntt::add_mod(u, v, mod);
            a[base + j + len / 2] = ntt::sub_mod(u, v, mod);
        }
        pthread_barrier_wait(&ctx.barrier);
    }

    if (ctx.inverse) {
        const u32 inv_n = ntt::pow_mod(static_cast<u32>(plan.n % mod), mod - 2u, mod);
        const std::size_t begin = plan.n * args->id / ctx.thread_count;
        const std::size_t end = plan.n * (args->id + 1) / ctx.thread_count;
        for (std::size_t i = begin; i < end; ++i) a[i] = ntt::mul_mod(a[i], inv_n, mod);
    }
    return nullptr;
}

void transform_pthread(std::vector<u32>& a, bool inverse, const Plan& plan, int threads) {
    ntt::apply_bit_reverse(a, plan);
    PthreadContext context{&a, &plan, inverse, threads, {}};
    if (pthread_barrier_init(&context.barrier, nullptr, static_cast<unsigned>(threads)) != 0) {
        throw std::runtime_error("pthread_barrier_init failed");
    }
    std::vector<pthread_t> handles(threads);
    std::vector<WorkerArgs> args(threads);
    for (int id = 0; id < threads; ++id) {
        args[id] = WorkerArgs{&context, id};
        if (pthread_create(&handles[id], nullptr, pthread_worker, &args[id]) != 0) {
            std::terminate();
        }
    }
    for (pthread_t handle : handles) pthread_join(handle, nullptr);
    pthread_barrier_destroy(&context.barrier);
}

void transform_openmp(std::vector<u32>& a, bool inverse, const Plan& plan, int threads) {
#ifndef _OPENMP
    (void)a; (void)inverse; (void)plan; (void)threads;
    throw std::runtime_error("OpenMP backend was not enabled at compile time");
#else
    const auto& stages = inverse ? plan.inverse_twiddles : plan.forward_twiddles;
    const u32 mod = plan.modulus.value;
    omp_set_dynamic(0);
    omp_set_num_threads(threads);

#pragma omp parallel for schedule(static)
    for (std::int64_t i = 0; i < static_cast<std::int64_t>(plan.n); ++i) {
        const std::size_t reversed = plan.bit_reverse[static_cast<std::size_t>(i)];
        if (static_cast<std::size_t>(i) < reversed) std::swap(a[static_cast<std::size_t>(i)], a[reversed]);
    }

#pragma omp parallel
    {
        for (std::size_t stage = 0, len = 2; stage < stages.size(); ++stage, len <<= 1u) {
            const auto& twiddle = stages[stage];
#pragma omp for schedule(static)
            for (std::int64_t butterfly = 0; butterfly < static_cast<std::int64_t>(plan.n / 2); ++butterfly) {
                const std::size_t b = static_cast<std::size_t>(butterfly);
                const std::size_t j = b % (len / 2);
                const std::size_t base = (b / (len / 2)) * len;
                const u32 u = a[base + j];
                const u32 v = ntt::mul_mod(a[base + j + len / 2], twiddle[j], mod);
                a[base + j] = ntt::add_mod(u, v, mod);
                a[base + j + len / 2] = ntt::sub_mod(u, v, mod);
            }
        }
        if (inverse) {
            const u32 inv_n = ntt::pow_mod(static_cast<u32>(plan.n % mod), mod - 2u, mod);
#pragma omp for schedule(static)
            for (std::int64_t i = 0; i < static_cast<std::int64_t>(plan.n); ++i) {
                a[static_cast<std::size_t>(i)] = ntt::mul_mod(a[static_cast<std::size_t>(i)], inv_n, mod);
            }
        }
    }
#endif
}

void process_block_simd(
    std::vector<u32>& a,
    std::size_t base,
    std::size_t len,
    const std::vector<u32>& twiddle,
    u32 mod,
    std::size_t begin = 0,
    std::size_t end = 0) {
    if (end == 0) end = len / 2;
    std::size_t j = begin;
#ifdef __AVX2__
    const __m256i v_mod = _mm256_set1_epi32(static_cast<int>(mod));
    const __m256i v_mod_minus_one = _mm256_set1_epi32(static_cast<int>(mod - 1));
    const __m256i v_zero = _mm256_setzero_si256();
    alignas(32) u32 products[8];
    for (; j + 8 <= end; j += 8) {
        for (std::size_t lane = 0; lane < 8; ++lane) {
            products[lane] = ntt::mul_mod(a[base + j + lane + len / 2], twiddle[j + lane], mod);
        }
        const __m256i u = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(&a[base + j]));
        const __m256i v = _mm256_load_si256(reinterpret_cast<const __m256i*>(products));
        __m256i sum = _mm256_add_epi32(u, v);
        const __m256i sum_mask = _mm256_cmpgt_epi32(sum, v_mod_minus_one);
        sum = _mm256_sub_epi32(sum, _mm256_and_si256(sum_mask, v_mod));
        __m256i difference = _mm256_sub_epi32(u, v);
        const __m256i diff_mask = _mm256_cmpgt_epi32(v_zero, difference);
        difference = _mm256_add_epi32(difference, _mm256_and_si256(diff_mask, v_mod));
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(&a[base + j]), sum);
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(&a[base + j + len / 2]), difference);
    }
#endif
    for (; j < end; ++j) {
        const u32 u = a[base + j];
        const u32 v = ntt::mul_mod(a[base + j + len / 2], twiddle[j], mod);
        a[base + j] = ntt::add_mod(u, v, mod);
        a[base + j + len / 2] = ntt::sub_mod(u, v, mod);
    }
}

void transform_simd(std::vector<u32>& a, bool inverse, const Plan& plan) {
    ntt::apply_bit_reverse(a, plan);
    const auto& stages = inverse ? plan.inverse_twiddles : plan.forward_twiddles;
    const u32 mod = plan.modulus.value;
    std::size_t len = 2;
    for (const auto& twiddle : stages) {
        for (std::size_t base = 0; base < plan.n; base += len) {
            process_block_simd(a, base, len, twiddle, mod);
        }
        len <<= 1u;
    }
    if (inverse) {
        const u32 inv_n = ntt::pow_mod(static_cast<u32>(plan.n % mod), mod - 2u, mod);
        for (u32& x : a) x = ntt::mul_mod(x, inv_n, mod);
    }
}

void transform_hybrid(std::vector<u32>& a, bool inverse, const Plan& plan, int threads) {
#ifndef _OPENMP
    (void)a; (void)inverse; (void)plan; (void)threads;
    throw std::runtime_error("hybrid backend requires OpenMP at compile time");
#else
    omp_set_dynamic(0);
    omp_set_num_threads(threads);
    const auto& stages = inverse ? plan.inverse_twiddles : plan.forward_twiddles;
    const u32 mod = plan.modulus.value;
#pragma omp parallel for schedule(static)
    for (std::int64_t i = 0; i < static_cast<std::int64_t>(plan.n); ++i) {
        const std::size_t reversed = plan.bit_reverse[static_cast<std::size_t>(i)];
        if (static_cast<std::size_t>(i) < reversed) std::swap(a[static_cast<std::size_t>(i)], a[reversed]);
    }
#pragma omp parallel
    {
        for (std::size_t stage = 0, len = 2; stage < stages.size(); ++stage, len <<= 1u) {
            const auto& twiddle = stages[stage];
            const std::size_t half = len / 2;
            const std::size_t chunks_per_group = (half + 7) / 8;
            const std::int64_t tasks = static_cast<std::int64_t>((plan.n / len) * chunks_per_group);
#pragma omp for schedule(static)
            for (std::int64_t task = 0; task < tasks; ++task) {
                const std::size_t task_index = static_cast<std::size_t>(task);
                const std::size_t group = task_index / chunks_per_group;
                const std::size_t begin = (task_index % chunks_per_group) * 8;
                const std::size_t end = std::min(begin + 8, half);
                process_block_simd(a, group * len, len, twiddle, mod, begin, end);
            }
        }
        if (inverse) {
            const u32 inv_n = ntt::pow_mod(static_cast<u32>(plan.n % mod), mod - 2u, mod);
#pragma omp for simd schedule(static)
            for (std::int64_t i = 0; i < static_cast<std::int64_t>(plan.n); ++i) {
                a[static_cast<std::size_t>(i)] = ntt::mul_mod(a[static_cast<std::size_t>(i)], inv_n, mod);
            }
        }
    }
#endif
}

void transform(std::vector<u32>& a, bool inverse, const Plan& plan, Backend backend, int threads) {
    switch (backend) {
        case Backend::Serial: ntt::transform_serial(a, inverse, plan); break;
        case Backend::Pthread: transform_pthread(a, inverse, plan, threads); break;
        case Backend::OpenMP: transform_openmp(a, inverse, plan, threads); break;
        case Backend::Simd: transform_simd(a, inverse, plan); break;
        case Backend::Hybrid: transform_hybrid(a, inverse, plan, threads); break;
    }
}

std::vector<u32> convolution(
    const std::vector<u32>& lhs,
    const std::vector<u32>& rhs,
    Backend backend,
    int threads) {
    const std::size_t output_size = lhs.size() + rhs.size() - 1;
    const std::size_t n = ntt::ceil_power_of_two(output_size);
    Plan plan(n, ntt::kDefaultMod);
    std::vector<u32> a(n), b(n);
    std::copy(lhs.begin(), lhs.end(), a.begin());
    std::copy(rhs.begin(), rhs.end(), b.begin());
    transform(a, false, plan, backend, threads);
    transform(b, false, plan, backend, threads);
#ifdef _OPENMP
#pragma omp parallel for simd schedule(static) if(backend == Backend::OpenMP || backend == Backend::Hybrid)
#endif
    for (std::int64_t i = 0; i < static_cast<std::int64_t>(n); ++i) {
        a[static_cast<std::size_t>(i)] = ntt::mul_mod(
            a[static_cast<std::size_t>(i)], b[static_cast<std::size_t>(i)], ntt::kDefaultMod.value);
    }
    transform(a, true, plan, backend, threads);
    a.resize(output_size);
    return a;
}

void run_self_test(int threads) {
    const auto lhs = ntt::make_input(127, 1000, 20260712);
    const auto rhs = ntt::make_input(83, 1000, 20260713);
    const auto expected = ntt::convolution_naive(lhs, rhs, ntt::kDefaultMod.value);
    for (Backend backend : {Backend::Serial, Backend::Pthread, Backend::OpenMP, Backend::Simd, Backend::Hybrid}) {
#ifndef _OPENMP
        if (backend == Backend::OpenMP || backend == Backend::Hybrid) continue;
#endif
        const auto actual = convolution(lhs, rhs, backend, threads);
        if (actual != expected) throw std::runtime_error(std::string("self-test failed: ") + backend_name(backend));
        std::cout << "self-test," << backend_name(backend) << ",PASS\n";
    }
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const Options options = parse_options(argc, argv);
        if (options.self_test) {
            run_self_test(options.threads);
            return 0;
        }

        const std::size_t polynomial_size = options.transform_size / 2;
        const auto lhs = ntt::make_input(polynomial_size, 1000, 20260712);
        const auto rhs = ntt::make_input(polynomial_size, 1000, 20260713);
        const auto reference = convolution(lhs, rhs, Backend::Serial, 1);
        std::vector<double> samples;
        std::vector<u32> result;
        for (int repeat = 0; repeat < options.repeats; ++repeat) {
            const auto start = std::chrono::steady_clock::now();
            result = convolution(lhs, rhs, options.backend, options.threads);
            samples.push_back(ntt::milliseconds_since(start));
        }
        if (result != reference) throw std::runtime_error("benchmark result differs from serial reference");
        const double mean = std::accumulate(samples.begin(), samples.end(), 0.0) / samples.size();
        double variance = 0.0;
        for (double sample : samples) variance += (sample - mean) * (sample - mean);
        variance /= samples.size();

        std::cout << "backend,n,threads,repeats,mean_ms,stddev_ms,correct\n"
                  << backend_name(options.backend) << ',' << options.transform_size << ','
                  << options.threads << ',' << options.repeats << ',' << std::fixed << std::setprecision(3)
                  << mean << ',' << std::sqrt(variance) << ",true\n";
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return 1;
    }
    return 0;
}
