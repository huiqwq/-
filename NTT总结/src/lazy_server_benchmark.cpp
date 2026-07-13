#include "ntt_common.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace {

using ntt::Plan;
using ntt::u32;
using ntt::u64;
using ntt::u128;

constexpr u32 kMod = 998244353u;
constexpr u64 kBarrettMu = static_cast<u64>((static_cast<u128>(1) << 64u) / kMod);

enum class Variant { BarrettSerial, LazySerial, BarrettOpenMP, LazyOpenMP };

struct Options {
    std::size_t n = 1u << 17u;
    int threads = 4;
    int repeats = 10;
    int warmups = 2;
    Variant variant = Variant::LazyOpenMP;
    bool self_test = false;
};

const char* variant_name(Variant variant) {
    switch (variant) {
        case Variant::BarrettSerial: return "barrett-serial";
        case Variant::LazySerial: return "barrett-lazy-serial";
        case Variant::BarrettOpenMP: return "barrett-openmp";
        case Variant::LazyOpenMP: return "barrett-lazy-openmp";
    }
    return "unknown";
}

Variant parse_variant(std::string_view value) {
    if (value == "barrett-serial") return Variant::BarrettSerial;
    if (value == "barrett-lazy-serial") return Variant::LazySerial;
    if (value == "barrett-openmp") return Variant::BarrettOpenMP;
    if (value == "barrett-lazy-openmp") return Variant::LazyOpenMP;
    throw std::invalid_argument("invalid --variant");
}

Options parse_options(int argc, char** argv) {
    Options options;
    for (int i = 1; i < argc; ++i) {
        const std::string_view arg(argv[i]);
        auto value = [&]() -> const char* {
            if (++i == argc) throw std::invalid_argument("missing command-line value");
            return argv[i];
        };
        if (arg == "--variant") options.variant = parse_variant(value());
        else if (arg == "--size") options.n = std::stoull(value());
        else if (arg == "--threads") options.threads = std::stoi(value());
        else if (arg == "--repeats") options.repeats = std::stoi(value());
        else if (arg == "--warmups") options.warmups = std::stoi(value());
        else if (arg == "--self-test") options.self_test = true;
        else if (arg == "--help") {
            std::cout
                << "Usage: ntt_lazy_server --variant barrett-serial|barrett-lazy-serial|"
                   "barrett-openmp|barrett-lazy-openmp [--size N] [--threads T] "
                   "[--warmups W] [--repeats R] [--self-test]\n";
            std::exit(0);
        } else {
            throw std::invalid_argument("unknown argument: " + std::string(arg));
        }
    }
    if (!ntt::is_power_of_two(options.n) || options.n < 2 || options.n > (1u << 23u)) {
        throw std::invalid_argument("--size must be a power of two in [2, 2^23]");
    }
    if (options.threads < 1 || options.repeats < 1 || options.warmups < 0) {
        throw std::invalid_argument("invalid threads/repeats/warmups");
    }
    return options;
}

inline u32 barrett_mul(u32 left, u32 right) {
    const u64 product = static_cast<u64>(left) * right;
    const u64 quotient = static_cast<u64>((static_cast<u128>(product) * kBarrettMu) >> 64u);
    u64 remainder = product - quotient * kMod;
    if (remainder >= kMod) remainder -= kMod;
    return static_cast<u32>(remainder);
}

inline u32 normalize(u32 value) {
    constexpr u32 two_mod = 2u * kMod;
    if (value >= two_mod) value -= two_mod;
    if (value >= kMod) value -= kMod;
    return value;
}

void transform_barrett(std::vector<u32>& data, bool inverse, const Plan& plan,
                       bool lazy, bool parallel, int threads) {
    const std::int64_t count = static_cast<std::int64_t>(plan.n);
#ifdef _OPENMP
    if (parallel) {
        omp_set_dynamic(0);
        omp_set_num_threads(threads);
    }
#else
    if (parallel) throw std::runtime_error("OpenMP variant requires -fopenmp");
    (void)threads;
#endif

    auto reverse_one = [&](std::int64_t i) {
        const std::size_t index = static_cast<std::size_t>(i);
        const std::size_t reversed = plan.bit_reverse[index];
        if (index < reversed) std::swap(data[index], data[reversed]);
    };
#ifdef _OPENMP
    if (parallel) {
#pragma omp parallel for schedule(static)
        for (std::int64_t i = 0; i < count; ++i) reverse_one(i);
    } else
#endif
    {
        for (std::int64_t i = 0; i < count; ++i) reverse_one(i);
    }

    const auto& stages = inverse ? plan.inverse_twiddles : plan.forward_twiddles;
    std::size_t length = 2;
    for (const auto& twiddles : stages) {
        const std::size_t half = length / 2;
        auto butterfly_one = [&](std::int64_t butterfly) {
            const std::size_t b = static_cast<std::size_t>(butterfly);
            const std::size_t j = b % half;
            const std::size_t first = (b / half) * length + j;
            const std::size_t second = first + half;
            u32 left = data[first];
            const u32 right = barrett_mul(data[second], twiddles[j]);
            if (lazy) {
                constexpr u32 two_mod = 2u * kMod;
                if (left >= two_mod) left -= two_mod;
                data[first] = left + right;
                data[second] = left + two_mod - right;
            } else {
                const u32 sum = left + right;
                data[first] = sum >= kMod ? sum - kMod : sum;
                data[second] = left >= right ? left - right : left + kMod - right;
            }
        };
#ifdef _OPENMP
        if (parallel) {
#pragma omp parallel for schedule(static)
            for (std::int64_t butterfly = 0; butterfly < count / 2; ++butterfly) butterfly_one(butterfly);
        } else
#endif
        {
            for (std::int64_t butterfly = 0; butterfly < count / 2; ++butterfly) butterfly_one(butterfly);
        }
        length <<= 1u;
    }

    if (inverse) {
        const u32 inverse_n = ntt::pow_mod(static_cast<u32>(plan.n), kMod - 2u, kMod);
        auto scale_one = [&](std::int64_t i) {
            const std::size_t index = static_cast<std::size_t>(i);
            data[index] = barrett_mul(lazy ? normalize(data[index]) : data[index], inverse_n);
        };
#ifdef _OPENMP
        if (parallel) {
#pragma omp parallel for schedule(static)
            for (std::int64_t i = 0; i < count; ++i) scale_one(i);
        } else
#endif
        {
            for (std::int64_t i = 0; i < count; ++i) scale_one(i);
        }
    }
}

void pointwise(std::vector<u32>& left, const std::vector<u32>& right, bool parallel) {
    const std::int64_t count = static_cast<std::int64_t>(left.size());
    auto multiply_one = [&](std::int64_t i) {
        const std::size_t index = static_cast<std::size_t>(i);
        left[index] = barrett_mul(left[index], right[index]);
    };
#ifdef _OPENMP
    if (parallel) {
#pragma omp parallel for schedule(static)
        for (std::int64_t i = 0; i < count; ++i) multiply_one(i);
    } else
#endif
    {
        for (std::int64_t i = 0; i < count; ++i) multiply_one(i);
    }
}

std::vector<u32> convolution(const std::vector<u32>& lhs, const std::vector<u32>& rhs,
                             const Plan& plan, Variant variant) {
    const bool lazy = variant == Variant::LazySerial || variant == Variant::LazyOpenMP;
    const bool parallel = variant == Variant::BarrettOpenMP || variant == Variant::LazyOpenMP;
    int threads = 1;
#ifdef _OPENMP
    threads = omp_get_max_threads();
#endif
    std::vector<u32> left(plan.n, 0), right(plan.n, 0);
    std::copy(lhs.begin(), lhs.end(), left.begin());
    std::copy(rhs.begin(), rhs.end(), right.begin());
    transform_barrett(left, false, plan, lazy, parallel, threads);
    transform_barrett(right, false, plan, lazy, parallel, threads);
    pointwise(left, right, parallel);
    transform_barrett(left, true, plan, lazy, parallel, threads);
    left.resize(lhs.size() + rhs.size() - 1);
    return left;
}

struct Sample { double core_ms; double end_ms; };

Sample run_once(const std::vector<u32>& lhs, const std::vector<u32>& rhs,
                const Plan& plan, Variant variant, std::vector<u32>& output) {
    const auto end_start = std::chrono::steady_clock::now();
    std::vector<u32> left(plan.n, 0), right(plan.n, 0);
    std::copy(lhs.begin(), lhs.end(), left.begin());
    std::copy(rhs.begin(), rhs.end(), right.begin());
    const auto core_start = std::chrono::steady_clock::now();

    const bool lazy = variant == Variant::LazySerial || variant == Variant::LazyOpenMP;
    const bool parallel = variant == Variant::BarrettOpenMP || variant == Variant::LazyOpenMP;
    int threads = 1;
#ifdef _OPENMP
    threads = omp_get_max_threads();
#endif
    transform_barrett(left, false, plan, lazy, parallel, threads);
    transform_barrett(right, false, plan, lazy, parallel, threads);
    pointwise(left, right, parallel);
    transform_barrett(left, true, plan, lazy, parallel, threads);
    const double core_ms = ntt::milliseconds_since(core_start);
    left.resize(lhs.size() + rhs.size() - 1);
    output = std::move(left);
    return {core_ms, ntt::milliseconds_since(end_start)};
}

double mean(const std::vector<double>& values) {
    return std::accumulate(values.begin(), values.end(), 0.0) / values.size();
}

double stddev(const std::vector<double>& values, double average) {
    double sum = 0.0;
    for (double value : values) sum += (value - average) * (value - average);
    return std::sqrt(sum / values.size());
}

void self_test(int threads) {
#ifdef _OPENMP
    omp_set_num_threads(threads);
#else
    (void)threads;
#endif
    const auto lhs = ntt::make_input(127, 1000, 20260712);
    const auto rhs = ntt::make_input(83, 1000, 20260713);
    const Plan plan(ntt::ceil_power_of_two(lhs.size() + rhs.size() - 1), ntt::kDefaultMod);
    const auto expected = ntt::convolution_naive(lhs, rhs, kMod);
    for (Variant variant : {Variant::BarrettSerial, Variant::LazySerial,
                            Variant::BarrettOpenMP, Variant::LazyOpenMP}) {
#ifndef _OPENMP
        if (variant == Variant::BarrettOpenMP || variant == Variant::LazyOpenMP) continue;
#endif
        if (convolution(lhs, rhs, plan, variant) != expected) {
            throw std::runtime_error(std::string("self-test failed: ") + variant_name(variant));
        }
        std::cout << "self-test," << variant_name(variant) << ",PASS\n";
    }
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const Options options = parse_options(argc, argv);
#ifdef _OPENMP
        omp_set_dynamic(0);
        omp_set_num_threads(options.threads);
#endif
        if (options.self_test) {
            self_test(options.threads);
            return 0;
        }
        const std::size_t input_size = options.n / 2;
        const auto lhs = ntt::make_input(input_size, 1000, 20260712);
        const auto rhs = ntt::make_input(input_size, 1000, 20260713);
        const Plan plan(options.n, ntt::kDefaultMod);  // deliberately outside timing
        const auto expected = ntt::convolution_serial(lhs, rhs, ntt::kDefaultMod);
        std::vector<u32> output;
        for (int i = 0; i < options.warmups; ++i) run_once(lhs, rhs, plan, options.variant, output);

        std::vector<double> core_samples, end_samples;
        for (int i = 0; i < options.repeats; ++i) {
            const Sample sample = run_once(lhs, rhs, plan, options.variant, output);
            core_samples.push_back(sample.core_ms);
            end_samples.push_back(sample.end_ms);
        }
        const bool correct = output == expected;
        const double core_mean = mean(core_samples);
        const double end_mean = mean(end_samples);
        std::cout << "platform,variant,n,threads_or_block,repeats,core_ms,end_to_end_ms,stddev_ms,correct\n"
                  << "server," << variant_name(options.variant) << ',' << options.n << ','
                  << options.threads << ',' << options.repeats << ',' << std::fixed
                  << std::setprecision(4) << core_mean << ',' << end_mean << ','
                  << stddev(core_samples, core_mean) << ',' << (correct ? "true" : "false") << '\n';
        if (!correct) return 2;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return 1;
    }
    return 0;
}
