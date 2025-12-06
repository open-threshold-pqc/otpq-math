#pragma once

/**
 * @file prime_utils.h
 * @brief Compile-time and runtime primality utilities for small and medium integers.
 *
 * These helpers are lightweight, constexpr-capable, and suitable for use in
 * template-based algebraic structures such as PrimeField or NTT parameter checks.
 *
 * The primality test implemented here is:
 *   - deterministic for 32-bit integers,
 *   - efficient for typical cryptographic moduli used in ℤ_p arithmetic,
 *   - constexpr-evaluatable for compile-time modulus validation.
 */

#include <concepts>
#include <cstdint>

namespace otpq::utils::math {
    /**
     * @brief Deterministic primality test for integral types.
     *
     * Works constexpr and is suitable for:
     *   - validating modulus `p` for finite fields 𝔽ₚ,
     *   - small-range compile-time checks,
     *   - runtime checks for typical moduli (≤ 2³¹).
     *
     * @tparam T Integral type (signed or unsigned).
     * @param n The number to test for primality.
     * @return true if n is prime, false otherwise.
     *
     * @note Uses the standard 6k ± 1 optimization.
     */
    template<std::integral T>
    constexpr bool is_prime(T n) noexcept {
        if (n <= 1) return false;
        if (n <= 3) return true;
        if (n % 2 == 0 || n % 3 == 0) return false;

        for (T i = 5; i * i <= n; i += 6) {
            if (n % i == 0 || n % (i + 2) == 0)
                return false;
        }
        return true;
    }
}
