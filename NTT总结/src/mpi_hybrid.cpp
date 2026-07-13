#include "ntt_common.hpp"

#include <mpi.h>

#include <array>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <string_view>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace {

using ntt::Plan;
using ntt::u32;
using ntt::u64;
using ntt::u128;

constexpr std::array<ntt::Modulus, 4> kCrtModuli{{
    {998244353u, 3u},
    {1004535809u, 3u},
    {469762049u, 3u},
    {1224736769u, 3u},
}};

struct Options {
    std::size_t transform_size = 1u << 17u;
    int repeats = 5;
    int threads = 2;
    u64 target_modulus = 1337006139375617ull;
    u32 coefficient_limit = 1000;
};

Options parse_options(int argc, char** argv) {
    Options options;
    for (int i = 1; i < argc; ++i) {
        const std::string_view arg(argv[i]);
        auto value = [&]() -> const char* {
            if (++i == argc) throw std::invalid_argument("missing command-line value");
            return argv[i];
        };
        if (arg == "--size") options.transform_size = std::stoull(value());
        else if (arg == "--repeats") options.repeats = std::stoi(value());
        else if (arg == "--threads") options.threads = std::stoi(value());
        else if (arg == "--target-mod") options.target_modulus = std::stoull(value());
        else if (arg == "--coefficient-limit") options.coefficient_limit = static_cast<u32>(std::stoul(value()));
        else if (arg == "--help") {
            std::cout << "Usage: mpiexec -np 4 ntt_mpi [--size N] [--threads T] "
                         "[--repeats R] [--target-mod P] [--coefficient-limit C]\n";
            std::exit(0);
        } else {
            throw std::invalid_argument("unknown argument: " + std::string(arg));
        }
    }
    if (!ntt::is_power_of_two(options.transform_size) || options.transform_size < 2) {
        throw std::invalid_argument("--size must be a power of two >= 2");
    }
    if (options.transform_size > (1u << 21u)) {
        throw std::invalid_argument("the selected CRT moduli support at most 2^21 points");
    }
    if (options.repeats < 1 || options.threads < 1 || options.target_modulus < 2 || options.coefficient_limit < 1) {
        throw std::invalid_argument("invalid non-positive option");
    }
    return options;
}

void transform_local(std::vector<u32>& a, bool inverse, const Plan& plan, int threads) {
    const auto& stages = inverse ? plan.inverse_twiddles : plan.forward_twiddles;
    const u32 mod = plan.modulus.value;
#ifdef _OPENMP
    omp_set_dynamic(0);
    omp_set_num_threads(threads);
#pragma omp parallel for schedule(static)
#endif
    for (std::int64_t i = 0; i < static_cast<std::int64_t>(plan.n); ++i) {
        const std::size_t reversed = plan.bit_reverse[static_cast<std::size_t>(i)];
        if (static_cast<std::size_t>(i) < reversed) std::swap(a[static_cast<std::size_t>(i)], a[reversed]);
    }

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        for (std::size_t stage = 0, len = 2; stage < stages.size(); ++stage, len <<= 1u) {
            const auto& twiddle = stages[stage];
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
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
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
            for (std::int64_t i = 0; i < static_cast<std::int64_t>(plan.n); ++i) {
                a[static_cast<std::size_t>(i)] = ntt::mul_mod(a[static_cast<std::size_t>(i)], inv_n, mod);
            }
        }
    }
}

std::vector<u32> local_convolution(
    const std::vector<u32>& lhs,
    const std::vector<u32>& rhs,
    std::size_t transform_size,
    ntt::Modulus modulus,
    int threads) {
    Plan plan(transform_size, modulus);
    std::vector<u32> a(transform_size), b(transform_size);
    for (std::size_t i = 0; i < lhs.size(); ++i) {
        a[i] = lhs[i] % modulus.value;
        b[i] = rhs[i] % modulus.value;
    }
    transform_local(a, false, plan, threads);
    transform_local(b, false, plan, threads);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (std::int64_t i = 0; i < static_cast<std::int64_t>(transform_size); ++i) {
        a[static_cast<std::size_t>(i)] = ntt::mul_mod(
            a[static_cast<std::size_t>(i)], b[static_cast<std::size_t>(i)], modulus.value);
    }
    transform_local(a, true, plan, threads);
    a.resize(lhs.size() + rhs.size() - 1);
    return a;
}

u64 garner_coefficient(const u32* residues, u64 target_modulus) {
    u128 x = residues[0];
    u128 product = kCrtModuli[0].value;
    for (std::size_t i = 1; i < kCrtModuli.size(); ++i) {
        const u32 prime = kCrtModuli[i].value;
        const u32 current = static_cast<u32>(x % prime);
        const u32 delta = ntt::sub_mod(residues[i], current, prime);
        const u32 product_mod = static_cast<u32>(product % prime);
        const u32 inverse = ntt::pow_mod(product_mod, prime - 2u, prime);
        const u32 digit = ntt::mul_mod(delta, inverse, prime);
        x += product * digit;
        product *= prime;
    }
    return static_cast<u64>(x % target_modulus);
}

u64 exact_coefficient(
    const std::vector<u32>& lhs,
    const std::vector<u32>& rhs,
    std::size_t index,
    u64 target_modulus) {
    const std::size_t begin = index >= rhs.size() - 1 ? index - (rhs.size() - 1) : 0;
    const std::size_t end = std::min(index, lhs.size() - 1);
    u128 sum = 0;
    for (std::size_t i = begin; i <= end; ++i) sum += static_cast<u128>(lhs[i]) * rhs[index - i];
    return static_cast<u64>(sum % target_modulus);
}

bool crt_range_is_safe(const Options& options) {
    long double crt_product = 1.0L;
    for (const auto modulus : kCrtModuli) crt_product *= modulus.value;
    const long double coefficient_bound = static_cast<long double>(options.transform_size / 2)
        * (options.coefficient_limit - 1.0L) * (options.coefficient_limit - 1.0L);
    return coefficient_bound < crt_product;
}

}  // namespace

int main(int argc, char** argv) {
    int provided = 0;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
    int rank = 0, world = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world);

    try {
        if (world != 4) throw std::runtime_error("this experiment requires exactly four MPI ranks");
        if (provided < MPI_THREAD_FUNNELED) throw std::runtime_error("MPI implementation lacks FUNNELED thread support");
        const Options options = parse_options(argc, argv);
        if (!crt_range_is_safe(options)) throw std::runtime_error("CRT product does not cover the coefficient bound");

        const std::size_t polynomial_size = options.transform_size / 2;
        const std::size_t output_size = polynomial_size * 2 - 1;
        if (output_size > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
            throw std::runtime_error("MPI count exceeds INT_MAX");
        }
        std::vector<u32> lhs(polynomial_size), rhs(polynomial_size);
        if (rank == 0) {
            lhs = ntt::make_input(polynomial_size, options.coefficient_limit, 20260712);
            rhs = ntt::make_input(polynomial_size, options.coefficient_limit, 20260713);
        }
        MPI_Bcast(lhs.data(), static_cast<int>(lhs.size()), MPI_UINT32_T, 0, MPI_COMM_WORLD);
        MPI_Bcast(rhs.data(), static_cast<int>(rhs.size()), MPI_UINT32_T, 0, MPI_COMM_WORLD);

        std::vector<double> elapsed_samples;
        bool correct = true;
        for (int repeat = 0; repeat < options.repeats; ++repeat) {
            MPI_Barrier(MPI_COMM_WORLD);
            const double begin = MPI_Wtime();
            const auto local = local_convolution(lhs, rhs, options.transform_size, kCrtModuli[rank], options.threads);
            std::vector<u32> gathered;
            if (rank == 0) gathered.resize(output_size * 4);
            MPI_Gather(
                local.data(), static_cast<int>(output_size), MPI_UINT32_T,
                rank == 0 ? gathered.data() : nullptr, static_cast<int>(output_size), MPI_UINT32_T,
                0, MPI_COMM_WORLD);

            if (rank == 0) {
                std::vector<u64> combined(output_size);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(options.threads)
#endif
                for (std::int64_t i = 0; i < static_cast<std::int64_t>(output_size); ++i) {
                    const std::size_t index = static_cast<std::size_t>(i);
                    const std::array<u32, 4> residues{{
                        gathered[index],
                        gathered[output_size + index],
                        gathered[2 * output_size + index],
                        gathered[3 * output_size + index],
                    }};
                    combined[index] = garner_coefficient(residues.data(), options.target_modulus);
                }
                for (std::size_t sample = 0; sample < 17; ++sample) {
                    const std::size_t index = sample * (output_size - 1) / 16;
                    if (combined[index] != exact_coefficient(lhs, rhs, index, options.target_modulus)) correct = false;
                }
            }
            const double local_elapsed = MPI_Wtime() - begin;
            double max_elapsed = 0.0;
            MPI_Reduce(&local_elapsed, &max_elapsed, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
            if (rank == 0) elapsed_samples.push_back(max_elapsed * 1000.0);
        }

        if (rank == 0) {
            const double mean = std::accumulate(elapsed_samples.begin(), elapsed_samples.end(), 0.0)
                / elapsed_samples.size();
            double variance = 0.0;
            for (double value : elapsed_samples) variance += (value - mean) * (value - mean);
            variance /= elapsed_samples.size();
            std::cout << "backend,n,ranks,threads_per_rank,repeats,mean_ms,stddev_ms,correct\n"
                      << "mpi+openmp+simd," << options.transform_size << ",4," << options.threads << ','
                      << options.repeats << ',' << std::fixed << std::setprecision(3) << mean << ','
                      << std::sqrt(variance) << ',' << (correct ? "true" : "false") << '\n';
            if (!correct) throw std::runtime_error("sampled CRT correctness check failed");
        }
    } catch (const std::exception& error) {
        std::cerr << "rank " << rank << " error: " << error.what() << '\n';
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    MPI_Finalize();
    return 0;
}
