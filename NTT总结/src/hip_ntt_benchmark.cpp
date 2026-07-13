#include "ntt_common.hpp"

#include <hip/hip_runtime.h>

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

namespace {

using ntt::u32;
using ntt::u64;

constexpr u32 kMod = 998244353u;
constexpr u32 kRoot = 3u;
constexpr u64 kBarrettMu = 18446744073709551615ull / kMod;
constexpr u32 kMontgomeryNInv = 998244351u;
constexpr u32 kMontgomeryR = 301989884u;
constexpr u32 kMontgomeryR2 = 932051910u;

enum Mode : int { Runtime = 0, Barrett = 1, Montgomery = 2, Lazy = 3 };
enum class Variant { Runtime, Barrett, Montgomery, Tiled, Lazy };

#define HIP_CHECK(call)                                                                        \
    do {                                                                                       \
        const hipError_t error__ = (call);                                                     \
        if (error__ != hipSuccess) {                                                           \
            throw std::runtime_error(std::string(#call) + ": " + hipGetErrorString(error__)); \
        }                                                                                      \
    } while (false)

struct Options {
    std::size_t n = 1u << 20u;
    int block_size = 128;
    int repeats = 10;
    int warmups = 2;
    Variant variant = Variant::Barrett;
    bool self_test = false;
};

const char* variant_name(Variant variant) {
    switch (variant) {
        case Variant::Runtime: return "runtime";
        case Variant::Barrett: return "barrett";
        case Variant::Montgomery: return "montgomery";
        case Variant::Tiled: return "montgomery-tiled";
        case Variant::Lazy: return "barrett-lazy";
    }
    return "unknown";
}

Variant parse_variant(std::string_view value) {
    if (value == "runtime") return Variant::Runtime;
    if (value == "barrett") return Variant::Barrett;
    if (value == "montgomery") return Variant::Montgomery;
    if (value == "tiled" || value == "montgomery-tiled") return Variant::Tiled;
    if (value == "lazy" || value == "barrett-lazy") return Variant::Lazy;
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
        else if (arg == "--block-size") options.block_size = std::stoi(value());
        else if (arg == "--repeats") options.repeats = std::stoi(value());
        else if (arg == "--warmups") options.warmups = std::stoi(value());
        else if (arg == "--self-test") options.self_test = true;
        else if (arg == "--help") {
            std::cout << "Usage: hip_ntt_benchmark --variant runtime|barrett|montgomery|tiled|lazy "
                         "[--size N] [--block-size B] [--warmups W] [--repeats R] [--self-test]\n";
            std::exit(0);
        } else {
            throw std::invalid_argument("unknown argument: " + std::string(arg));
        }
    }
    if (!ntt::is_power_of_two(options.n) || options.n < 2 || options.n > (1u << 23u)) {
        throw std::invalid_argument("--size must be a power of two in [2, 2^23]");
    }
    if (options.block_size < 32 || options.block_size > 1024 ||
        options.repeats < 1 || options.warmups < 0) {
        throw std::invalid_argument("invalid block-size/repeats/warmups");
    }
    return options;
}

__host__ __device__ __forceinline__ u32 montgomery_reduce(u64 value) {
    const u32 multiplier = static_cast<u32>(value) * kMontgomeryNInv;
    u64 reduced = (value + static_cast<u64>(multiplier) * kMod) >> 32u;
    if (reduced >= kMod) reduced -= kMod;
    return static_cast<u32>(reduced);
}

__host__ __device__ __forceinline__ u32 montgomery_mul(u32 left, u32 right) {
    return montgomery_reduce(static_cast<u64>(left) * right);
}

__device__ __forceinline__ u32 barrett_mul(u32 left, u32 right) {
    const u64 product = static_cast<u64>(left) * right;
    const u64 quotient = __umul64hi(static_cast<unsigned long long>(product),
                                    static_cast<unsigned long long>(kBarrettMu));
    u64 remainder = product - quotient * static_cast<u64>(kMod);
    if (remainder >= kMod) remainder -= kMod;
    return static_cast<u32>(remainder);
}

template <int mode>
__device__ __forceinline__ u32 multiply(u32 left, u32 right, u32 runtime_mod) {
    if constexpr (mode == Runtime) {
        return static_cast<u32>(static_cast<u64>(left) * right % runtime_mod);
    } else if constexpr (mode == Montgomery) {
        return montgomery_mul(left, right);
    } else {
        return barrett_mul(left, right);
    }
}

__global__ void bit_reverse_kernel(u32* data, const u32* reversed, std::size_t n) {
    const std::size_t index = static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
    if (index >= n) return;
    const std::size_t other = reversed[index];
    if (index < other) {
        const u32 temporary = data[index];
        data[index] = data[other];
        data[other] = temporary;
    }
}

__global__ void to_montgomery_kernel(u32* data, std::size_t n) {
    const std::size_t index = static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
    if (index < n) data[index] = montgomery_mul(data[index], kMontgomeryR2);
}

template <int mode>
__global__ void butterfly_kernel(u32* data, const u32* twiddles,
                                 std::size_t n, std::size_t length, u32 runtime_mod) {
    const std::size_t butterfly = static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
    if (butterfly >= n / 2) return;
    const std::size_t half = length / 2;
    const std::size_t group = butterfly / half;
    const std::size_t j = butterfly % half;
    const std::size_t first = group * length + j;
    const std::size_t second = first + half;
    u32 left = data[first];
    const u32 right = multiply<mode>(data[second], twiddles[j * (n / length)], runtime_mod);

    if constexpr (mode == Lazy) {
        constexpr u32 two_mod = 2u * kMod;
        if (left >= two_mod) left -= two_mod;
        data[first] = left + right;
        data[second] = left + two_mod - right;
    } else {
        const u32 sum = left + right;
        data[first] = sum >= kMod ? sum - kMod : sum;
        data[second] = left >= right ? left - right : left + kMod - right;
    }
}

template <int mode>
__global__ void pointwise_kernel(u32* left, const u32* right, std::size_t n, u32 runtime_mod) {
    const std::size_t index = static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
    if (index < n) left[index] = multiply<mode>(left[index], right[index], runtime_mod);
}

template <int mode>
__global__ void scale_kernel(u32* data, std::size_t n, u32 inverse_n, u32 runtime_mod) {
    const std::size_t index = static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
    if (index >= n) return;
    u32 value = data[index];
    if constexpr (mode == Lazy) {
        constexpr u32 two_mod = 2u * kMod;
        if (value >= two_mod) value -= two_mod;
        if (value >= kMod) value -= kMod;
    }
    value = multiply<mode>(value, inverse_n, runtime_mod);
    if constexpr (mode == Montgomery) value = montgomery_reduce(value);
    data[index] = value;
}

// DIT stages length=2..256 are local after the global bit-reversal permutation.
__global__ void tiled_first_eight_kernel(u32* data, const u32* twiddles, std::size_t n) {
    __shared__ u32 tile[256];
    const unsigned lane = threadIdx.x;
    const std::size_t base = static_cast<std::size_t>(blockIdx.x) * 256u;
    if (lane < 128u) {
        tile[lane] = data[base + lane];
        tile[lane + 128u] = data[base + lane + 128u];
    }
    __syncthreads();

    for (unsigned length = 2; length <= 256; length <<= 1u) {
        if (lane < 128u) {
            const unsigned half = length / 2u;
            const unsigned group = lane / half;
            const unsigned j = lane % half;
            const unsigned first = group * length + j;
            const unsigned second = first + half;
            const u32 left = tile[first];
            const u32 right = montgomery_mul(tile[second], twiddles[j * (n / length)]);
            const u32 sum = left + right;
            tile[first] = sum >= kMod ? sum - kMod : sum;
            tile[second] = left >= right ? left - right : left + kMod - right;
        }
        __syncthreads();
    }
    if (lane < 128u) {
        data[base + lane] = tile[lane];
        data[base + lane + 128u] = tile[lane + 128u];
    }
}

std::vector<u32> make_twiddles(std::size_t n, bool inverse, bool montgomery) {
    std::vector<u32> values(n / 2);
    u32 omega = ntt::pow_mod(kRoot, (kMod - 1u) / n, kMod);
    if (inverse) omega = ntt::pow_mod(omega, kMod - 2u, kMod);
    values[0] = 1;
    for (std::size_t i = 1; i < values.size(); ++i) values[i] = ntt::mul_mod(values[i - 1], omega, kMod);
    if (montgomery) {
        for (u32& value : values) value = static_cast<u32>(static_cast<u64>(value) * kMontgomeryR % kMod);
    }
    return values;
}

std::vector<u32> make_bit_reverse(std::size_t n) {
    std::vector<u32> table(n);
    unsigned bits = 0;
    while ((std::size_t{1} << bits) < n) ++bits;
    for (std::size_t i = 0; i < n; ++i) {
        std::size_t value = i, reversed = 0;
        for (unsigned bit = 0; bit < bits; ++bit) {
            reversed = (reversed << 1u) | (value & 1u);
            value >>= 1u;
        }
        table[i] = static_cast<u32>(reversed);
    }
    return table;
}

template <int mode>
void launch_stage(u32* data, const u32* twiddles, std::size_t n,
                  std::size_t length, int block_size, hipStream_t stream) {
    const std::size_t grid = (n / 2 + block_size - 1) / block_size;
    hipLaunchKernelGGL(HIP_KERNEL_NAME(butterfly_kernel<mode>), dim3(grid), dim3(block_size), 0,
                       stream, data, twiddles, n, length, kMod);
    HIP_CHECK(hipGetLastError());
}

void transform(u32* data, const u32* reversed, const u32* twiddles, std::size_t n,
               bool inverse, Variant variant, int block_size, hipStream_t stream) {
    const std::size_t element_grid = (n + block_size - 1) / block_size;
    hipLaunchKernelGGL(bit_reverse_kernel, dim3(element_grid), dim3(block_size), 0, stream,
                       data, reversed, n);
    HIP_CHECK(hipGetLastError());

    std::size_t first_standard_stage = 2;
    if (variant == Variant::Tiled && n >= 256) {
        hipLaunchKernelGGL(tiled_first_eight_kernel, dim3(n / 256), dim3(128), 0, stream,
                           data, twiddles, n);
        HIP_CHECK(hipGetLastError());
        first_standard_stage = 512;
    }
    for (std::size_t length = first_standard_stage; length <= n; length <<= 1u) {
        switch (variant) {
            case Variant::Runtime: launch_stage<Runtime>(data, twiddles, n, length, block_size, stream); break;
            case Variant::Barrett: launch_stage<Barrett>(data, twiddles, n, length, block_size, stream); break;
            case Variant::Lazy: launch_stage<Lazy>(data, twiddles, n, length, block_size, stream); break;
            case Variant::Montgomery:
            case Variant::Tiled: launch_stage<Montgomery>(data, twiddles, n, length, block_size, stream); break;
        }
    }

    if (inverse) {
        u32 inverse_n = ntt::pow_mod(static_cast<u32>(n), kMod - 2u, kMod);
        const int mode = (variant == Variant::Montgomery || variant == Variant::Tiled) ? Montgomery
                       : variant == Variant::Runtime ? Runtime : variant == Variant::Lazy ? Lazy : Barrett;
        if (mode == Montgomery) inverse_n = static_cast<u32>(static_cast<u64>(inverse_n) * kMontgomeryR % kMod);
        if (mode == Runtime) {
            hipLaunchKernelGGL(HIP_KERNEL_NAME(scale_kernel<Runtime>), dim3(element_grid), dim3(block_size),
                               0, stream, data, n, inverse_n, kMod);
        } else if (mode == Barrett) {
            hipLaunchKernelGGL(HIP_KERNEL_NAME(scale_kernel<Barrett>), dim3(element_grid), dim3(block_size),
                               0, stream, data, n, inverse_n, kMod);
        } else if (mode == Lazy) {
            hipLaunchKernelGGL(HIP_KERNEL_NAME(scale_kernel<Lazy>), dim3(element_grid), dim3(block_size),
                               0, stream, data, n, inverse_n, kMod);
        } else {
            hipLaunchKernelGGL(HIP_KERNEL_NAME(scale_kernel<Montgomery>), dim3(element_grid), dim3(block_size),
                               0, stream, data, n, inverse_n, kMod);
        }
        HIP_CHECK(hipGetLastError());
    }
}

struct DeviceBuffers {
    u32* left{};
    u32* right{};
    u32* reversed{};
    u32* forward{};
    u32* inverse{};
    hipStream_t stream{};
    hipEvent_t start{};
    hipEvent_t stop{};

    explicit DeviceBuffers(std::size_t n) {
        HIP_CHECK(hipMalloc(&left, n * sizeof(u32)));
        HIP_CHECK(hipMalloc(&right, n * sizeof(u32)));
        HIP_CHECK(hipMalloc(&reversed, n * sizeof(u32)));
        HIP_CHECK(hipMalloc(&forward, (n / 2) * sizeof(u32)));
        HIP_CHECK(hipMalloc(&inverse, (n / 2) * sizeof(u32)));
        HIP_CHECK(hipStreamCreate(&stream));
        HIP_CHECK(hipEventCreate(&start));
        HIP_CHECK(hipEventCreate(&stop));
    }
    ~DeviceBuffers() {
        if (start) hipEventDestroy(start);
        if (stop) hipEventDestroy(stop);
        if (stream) hipStreamDestroy(stream);
        if (left) hipFree(left);
        if (right) hipFree(right);
        if (reversed) hipFree(reversed);
        if (forward) hipFree(forward);
        if (inverse) hipFree(inverse);
    }
    DeviceBuffers(const DeviceBuffers&) = delete;
    DeviceBuffers& operator=(const DeviceBuffers&) = delete;
};

struct Sample { double core_ms; double end_ms; };

Sample run_once(const std::vector<u32>& lhs, const std::vector<u32>& rhs,
                std::vector<u32>& output, DeviceBuffers& device,
                std::size_t n, Variant variant, int block_size) {
    const auto end_start = std::chrono::steady_clock::now();
    HIP_CHECK(hipMemcpyAsync(device.left, lhs.data(), n * sizeof(u32), hipMemcpyHostToDevice, device.stream));
    HIP_CHECK(hipMemcpyAsync(device.right, rhs.data(), n * sizeof(u32), hipMemcpyHostToDevice, device.stream));
    HIP_CHECK(hipEventRecord(device.start, device.stream));
    const std::size_t grid = (n + block_size - 1) / block_size;
    if (variant == Variant::Montgomery || variant == Variant::Tiled) {
        hipLaunchKernelGGL(to_montgomery_kernel, dim3(grid), dim3(block_size), 0, device.stream,
                           device.left, n);
        hipLaunchKernelGGL(to_montgomery_kernel, dim3(grid), dim3(block_size), 0, device.stream,
                           device.right, n);
        HIP_CHECK(hipGetLastError());
    }
    transform(device.left, device.reversed, device.forward, n, false, variant, block_size, device.stream);
    transform(device.right, device.reversed, device.forward, n, false, variant, block_size, device.stream);

    if (variant == Variant::Runtime) {
        hipLaunchKernelGGL(HIP_KERNEL_NAME(pointwise_kernel<Runtime>), dim3(grid), dim3(block_size), 0,
                           device.stream, device.left, device.right, n, kMod);
    } else if (variant == Variant::Barrett) {
        hipLaunchKernelGGL(HIP_KERNEL_NAME(pointwise_kernel<Barrett>), dim3(grid), dim3(block_size), 0,
                           device.stream, device.left, device.right, n, kMod);
    } else if (variant == Variant::Lazy) {
        hipLaunchKernelGGL(HIP_KERNEL_NAME(pointwise_kernel<Lazy>), dim3(grid), dim3(block_size), 0,
                           device.stream, device.left, device.right, n, kMod);
    } else {
        hipLaunchKernelGGL(HIP_KERNEL_NAME(pointwise_kernel<Montgomery>), dim3(grid), dim3(block_size), 0,
                           device.stream, device.left, device.right, n, kMod);
    }
    HIP_CHECK(hipGetLastError());
    transform(device.left, device.reversed, device.inverse, n, true, variant, block_size, device.stream);
    HIP_CHECK(hipEventRecord(device.stop, device.stream));
    HIP_CHECK(hipMemcpyAsync(output.data(), device.left, n * sizeof(u32), hipMemcpyDeviceToHost, device.stream));
    HIP_CHECK(hipStreamSynchronize(device.stream));
    float core_ms = 0.0f;
    HIP_CHECK(hipEventElapsedTime(&core_ms, device.start, device.stop));
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

void upload_plan(DeviceBuffers& device, std::size_t n, Variant variant) {
    const bool montgomery = variant == Variant::Montgomery || variant == Variant::Tiled;
    const auto reversed = make_bit_reverse(n);
    const auto forward = make_twiddles(n, false, montgomery);
    const auto inverse = make_twiddles(n, true, montgomery);
    HIP_CHECK(hipMemcpy(device.reversed, reversed.data(), n * sizeof(u32), hipMemcpyHostToDevice));
    HIP_CHECK(hipMemcpy(device.forward, forward.data(), (n / 2) * sizeof(u32), hipMemcpyHostToDevice));
    HIP_CHECK(hipMemcpy(device.inverse, inverse.data(), (n / 2) * sizeof(u32), hipMemcpyHostToDevice));
}

bool benchmark(const Options& options, bool print_csv) {
    const std::size_t input_size = options.n / 2;
    std::vector<u32> lhs(options.n, 0), rhs(options.n, 0), output(options.n);
    const auto lhs_input = ntt::make_input(input_size, 1000, 20260712);
    const auto rhs_input = ntt::make_input(input_size, 1000, 20260713);
    std::copy(lhs_input.begin(), lhs_input.end(), lhs.begin());
    std::copy(rhs_input.begin(), rhs_input.end(), rhs.begin());
    const auto expected = ntt::convolution_serial(lhs_input, rhs_input, ntt::kDefaultMod);

    DeviceBuffers device(options.n);
    upload_plan(device, options.n, options.variant);  // plan and first transfer are outside timing
    for (int i = 0; i < options.warmups; ++i) {
        run_once(lhs, rhs, output, device, options.n, options.variant, options.block_size);
    }
    std::vector<double> core_samples, end_samples;
    for (int i = 0; i < options.repeats; ++i) {
        const Sample sample = run_once(lhs, rhs, output, device, options.n, options.variant, options.block_size);
        core_samples.push_back(sample.core_ms);
        end_samples.push_back(sample.end_ms);
    }
    const bool correct = std::equal(expected.begin(), expected.end(), output.begin());
    if (print_csv) {
        const double core_mean = mean(core_samples);
        std::cout << "platform,variant,n,threads_or_block,repeats,core_ms,end_to_end_ms,stddev_ms,correct\n"
                  << "hip," << variant_name(options.variant) << ',' << options.n << ','
                  << options.block_size << ',' << options.repeats << ',' << std::fixed
                  << std::setprecision(4) << core_mean << ',' << mean(end_samples) << ','
                  << stddev(core_samples, core_mean) << ',' << (correct ? "true" : "false") << '\n';
    }
    return correct;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const Options options = parse_options(argc, argv);
        if (options.self_test) {
            Options test = options;
            test.n = 1024;
            test.repeats = 1;
            test.warmups = 0;
            for (Variant variant : {Variant::Runtime, Variant::Barrett, Variant::Montgomery,
                                    Variant::Tiled, Variant::Lazy}) {
                test.variant = variant;
                if (!benchmark(test, false)) throw std::runtime_error(std::string("self-test failed: ") + variant_name(variant));
                std::cout << "self-test," << variant_name(variant) << ",PASS\n";
            }
            return 0;
        }
        if (!benchmark(options, true)) return 2;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return 1;
    }
    return 0;
}
