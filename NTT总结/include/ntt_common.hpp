#pragma once

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

namespace ntt {

using u32 = std::uint32_t;
using u64 = std::uint64_t;
using u128 = __uint128_t;

struct Modulus {
    u32 value;
    u32 primitive_root;
};

inline constexpr Modulus kDefaultMod{998244353u, 3u};

inline u32 add_mod(u32 a, u32 b, u32 mod) {
    const u32 sum = a + b;
    return sum >= mod ? sum - mod : sum;
}

inline u32 sub_mod(u32 a, u32 b, u32 mod) {
    return a >= b ? a - b : a + mod - b;
}

inline u32 mul_mod(u32 a, u32 b, u32 mod) {
    return static_cast<u32>(static_cast<u64>(a) * b % mod);
}

inline u32 pow_mod(u32 base, u64 exp, u32 mod) {
    u32 result = 1;
    while (exp != 0) {
        if (exp & 1u) result = mul_mod(result, base, mod);
        base = mul_mod(base, base, mod);
        exp >>= 1u;
    }
    return result;
}

inline bool is_power_of_two(std::size_t n) {
    return n != 0 && (n & (n - 1)) == 0;
}

inline std::size_t ceil_power_of_two(std::size_t n) {
    std::size_t result = 1;
    while (result < n) result <<= 1u;
    return result;
}

struct Plan {
    std::size_t n{};
    Modulus modulus{};
    std::vector<std::size_t> bit_reverse;
    std::vector<std::vector<u32>> forward_twiddles;
    std::vector<std::vector<u32>> inverse_twiddles;

    Plan(std::size_t size, Modulus mod) : n(size), modulus(mod) {
        if (!is_power_of_two(n)) throw std::invalid_argument("NTT length must be a power of two");
        if ((static_cast<u64>(modulus.value) - 1) % n != 0) {
            throw std::invalid_argument("NTT length must divide modulus - 1");
        }

        bit_reverse.resize(n);
        unsigned bits = 0;
        while ((std::size_t{1} << bits) < n) ++bits;
        for (std::size_t i = 0; i < n; ++i) {
            std::size_t x = i;
            std::size_t reversed = 0;
            for (unsigned b = 0; b < bits; ++b) {
                reversed = (reversed << 1u) | (x & 1u);
                x >>= 1u;
            }
            bit_reverse[i] = reversed;
        }

        const u32 inverse_root = pow_mod(modulus.primitive_root, modulus.value - 2u, modulus.value);
        for (std::size_t len = 2; len <= n; len <<= 1u) {
            forward_twiddles.push_back(make_stage_twiddles(len, modulus.primitive_root));
            inverse_twiddles.push_back(make_stage_twiddles(len, inverse_root));
        }
    }

private:
    std::vector<u32> make_stage_twiddles(std::size_t len, u32 root) const {
        const u32 stage_root = pow_mod(
            root, (static_cast<u64>(modulus.value) - 1) / len, modulus.value);
        std::vector<u32> table(len / 2);
        table[0] = 1;
        for (std::size_t j = 1; j < table.size(); ++j) {
            table[j] = mul_mod(table[j - 1], stage_root, modulus.value);
        }
        return table;
    }
};

inline void apply_bit_reverse(std::vector<u32>& a, const Plan& plan) {
    if (a.size() != plan.n) throw std::invalid_argument("plan/data size mismatch");
    for (std::size_t i = 0; i < plan.n; ++i) {
        if (i < plan.bit_reverse[i]) std::swap(a[i], a[plan.bit_reverse[i]]);
    }
}

inline void transform_serial(std::vector<u32>& a, bool inverse, const Plan& plan) {
    apply_bit_reverse(a, plan);
    const u32 mod = plan.modulus.value;
    const auto& stages = inverse ? plan.inverse_twiddles : plan.forward_twiddles;
    std::size_t len = 2;
    for (const auto& twiddle : stages) {
        for (std::size_t base = 0; base < plan.n; base += len) {
            for (std::size_t j = 0; j < len / 2; ++j) {
                const u32 u = a[base + j];
                const u32 v = mul_mod(a[base + j + len / 2], twiddle[j], mod);
                a[base + j] = add_mod(u, v, mod);
                a[base + j + len / 2] = sub_mod(u, v, mod);
            }
        }
        len <<= 1u;
    }
    if (inverse) {
        const u32 inv_n = pow_mod(static_cast<u32>(plan.n % mod), mod - 2u, mod);
        for (u32& x : a) x = mul_mod(x, inv_n, mod);
    }
}

inline std::vector<u32> convolution_serial(
    const std::vector<u32>& lhs,
    const std::vector<u32>& rhs,
    Modulus modulus = kDefaultMod) {
    if (lhs.empty() || rhs.empty()) return {};
    const std::size_t output_size = lhs.size() + rhs.size() - 1;
    const std::size_t n = ceil_power_of_two(output_size);
    Plan plan(n, modulus);
    std::vector<u32> a(n, 0), b(n, 0);
    std::copy(lhs.begin(), lhs.end(), a.begin());
    std::copy(rhs.begin(), rhs.end(), b.begin());
    transform_serial(a, false, plan);
    transform_serial(b, false, plan);
    for (std::size_t i = 0; i < n; ++i) a[i] = mul_mod(a[i], b[i], modulus.value);
    transform_serial(a, true, plan);
    a.resize(output_size);
    return a;
}

inline std::vector<u32> convolution_naive(
    const std::vector<u32>& lhs,
    const std::vector<u32>& rhs,
    u32 mod) {
    if (lhs.empty() || rhs.empty()) return {};
    std::vector<u32> result(lhs.size() + rhs.size() - 1, 0);
    for (std::size_t i = 0; i < lhs.size(); ++i) {
        for (std::size_t j = 0; j < rhs.size(); ++j) {
            const u32 product = mul_mod(lhs[i], rhs[j], mod);
            result[i + j] = add_mod(result[i + j], product, mod);
        }
    }
    return result;
}

inline std::vector<u32> make_input(std::size_t n, u32 limit, u64 seed) {
    if (limit == 0) throw std::invalid_argument("input limit must be positive");
    std::mt19937_64 generator(seed);
    std::uniform_int_distribution<u32> distribution(0, limit - 1);
    std::vector<u32> values(n);
    for (u32& value : values) value = distribution(generator);
    return values;
}

inline double milliseconds_since(std::chrono::steady_clock::time_point start) {
    return std::chrono::duration<double, std::milli>(
        std::chrono::steady_clock::now() - start).count();
}

}  // namespace ntt
