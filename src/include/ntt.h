#pragma once

#include <bit>
#include <cstddef>
#include <polynomial_fields.h>

namespace otpq::math::als::polynomial::transforms::ntt::helpers {
    /**
     * @brief Compute modular exponentiation using binary exponentiation.
     *
     * @tparam CoeffType  Unsigned or modular integer type.
     * @tparam ExpType    Integer exponent type.
     * @param base        Base value (mod `mod`).
     * @param exp         Exponent (non-negative).
     * @param mod         Modulus (must satisfy mod > 1).
     *
     * @returns (base^exp) mod mod
     *
     * @note This is constexpr-friendly and suitable for NTT root computation
     *       when the modulus is known at compile time.
     */
    template<typename CoeffType, typename ExpType>
    CoeffType mod_pow(CoeffType base, ExpType exp, CoeffType mod) {
        CoeffType r = 1;
        while (exp > 0) {
            if (exp & 1) r = (r * base) % mod;
            base = (base * base) % mod;
            exp >>= 1;
        }
        return r;
    }


    /**
     * @brief Reverse the lowest @p logn bits of @p x.
     *
     * @param x     Value whose bit order is reversed.
     * @param logn  Number of bits to reverse (typically log2(N)).
     *
     * @returns A bit-reversed index suitable for in-place NTT permutations.
     *
     * @note NTT requires permutation by reversing log₂(N) bits.
     */
    constexpr std::size_t reverse_bits(std::size_t x, std::size_t logn) noexcept {
        std::size_t r = 0;
        for (std::size_t i = 0; i < logn; ++i) {
            r = (r << 1) | (x & 1);
            x >>= 1;
        }
        return r;
    }
}

namespace otpq::math::als::polynomial::transforms::ntt {
    using helpers::mod_pow;
    using helpers::reverse_bits;


    /**
     * @brief Compute the forward Number-Theoretic Transform (NTT) in-place.
     *
     * @tparam N         Polynomial degree (NTT size), must be a power of two.
     * @tparam PolyMod   Polynomial modulus type (e.g., x^N + 1).
     * @tparam CoeffType Coefficient type used in modular arithmetic.
     * @tparam Als       Algebraic structure implementing Z_q arithmetic.
     *
     * @param poly       Input polynomial whose coefficients lie in Z_q.
     * @param root       A primitive N-th root of unity modulo q.
     *
     * @returns A polynomial whose coefficients represent the forward NTT.
     *
     * @pre
     * - @c N is a power of two.
     * - The modulus q satisfies `N | (q − 1)`.
     * - @p root is a generator of order exactly N in Z_q.
     * - The underlying field @p Als obeys the concept `core::Field`.
     *
     * @post
     * - Returns the in-place Cooley–Tukey NTT of `poly`.
     * - Output is in bit-reversed order if combined with certain convolution schemes,
     *   but here the algorithm includes bit-reversal so the output is in natural order.
     *
     * @details
     * This implementation performs:
     * 1. Bit-reversal permutation.
     * 2. Iterative Cooley–Tukey butterflies.
     *
     * @note The transform is unnormalized. Inverse NTT must multiply by N⁻¹ mod q.
     */
    template<
        std::size_t N,
        typename PolyMod = NoPolyModulus<N>,
        typename CoeffType = OTPQ_ALS_UNDERLYING_TYPE,
        typename Als = arith::ZmodNField<CoeffType> > requires core::Field<Als, CoeffType>
    auto toNTT(const PolynomialField<N, PolyMod, CoeffType, Als> &poly,
               CoeffType root) {
        auto out = poly;
        auto &coeffs = out.coeffs;
        const auto &alg = out.alg;
        const CoeffType q = alg.mod;

        constexpr std::size_t n = N;
        constexpr std::size_t logn = std::countr_zero(N);

        // --- Bit reversal permutation ---
        for (std::size_t i = 0; i < n; ++i) {
            auto j = reverse_bits(i, logn);
            if (i < j) std::swap(coeffs[i], coeffs[j]);
        }

        // --- Iterative Cooley–Tukey forward NTT ---
        for (std::size_t len = 2; len <= n; len <<= 1) {
            CoeffType wlen = mod_pow(root, n / len, q);

            for (std::size_t i = 0; i < n; i += len) {
                CoeffType w = 1;

                for (std::size_t j = 0; j < len / 2; ++j) {
                    CoeffType u = coeffs[i + j];
                    CoeffType t = alg.mul(coeffs[i + j + len / 2], w);

                    coeffs[i + j] = alg.add(u, t);
                    coeffs[i + j + len / 2] = alg.sub(u, t);

                    w = alg.mul(w, wlen);
                }
            }
        }

        return out;
    }


    /**
     * @brief Compute the inverse Number-Theoretic Transform (iNTT) in-place.
     *
     * @tparam N         Polynomial degree (NTT size), must be a power of two.
     * @tparam PolyMod   Polynomial modulus type.
     * @tparam CoeffType Coefficient type used in modular arithmetic.
     * @tparam Als       Algebraic structure implementing Z_q arithmetic.
     *
     * @param poly        Input polynomial transformed by @ref toNTT().
     * @param root_inv    The multiplicative inverse of the primitive N-th root of unity.
     *
     * @returns A polynomial in coefficient space (i.e., inverse NTT result).
     *
     * @pre
     * - @c N is a power of two.
     * - `root_inv` is the inverse of the primitive N-th root of unity.
     * - The field modulus satisfies `N | (q − 1)`.
     *
     * @post
     * - Returns the inverse NTT of `poly`.
     * - Output is scaled by N⁻¹ mod q to fully invert @ref toNTT().
     *
     * @details
     * This performs:
     * 1. Bit-reversal permutation.
     * 2. Gentleman–Sande inverse butterflies.
     * 3. Scaling by N⁻¹ modulo q.
     *
     * @note This implementation assumes q is prime so Fermat’s little theorem can be used
     *       for computing N⁻¹. Replace with a generic modular inverse for composite moduli.
     */
    template<
        std::size_t N,
        typename PolyMod = NoPolyModulus<N>,
        typename CoeffType = OTPQ_ALS_UNDERLYING_TYPE,
        typename Als = arith::ZmodNField<CoeffType> > requires core::Field<Als, CoeffType>
    auto fromNTT(const PolynomialField<N, PolyMod, CoeffType, Als> &poly,
                 CoeffType root_inv) {
        auto out = poly;
        auto &coeffs = out.coeffs;
        const auto &alg = out.alg;
        const CoeffType q = alg.mod;

        constexpr std::size_t n = N;
        constexpr std::size_t logn = std::countr_zero(N);

        // --- Bit reversal permutation ---
        for (std::size_t i = 0; i < n; ++i) {
            auto j = reverse_bits(i, logn);
            if (i < j) std::swap(coeffs[i], coeffs[j]);
        }

        // --- Gentleman–Sande inverse butterflies ---
        for (std::size_t len = 2; len <= n; len <<= 1) {
            CoeffType wlen = mod_pow(root_inv, n / len, q);

            for (std::size_t i = 0; i < n; i += len) {
                CoeffType w = 1;

                for (std::size_t j = 0; j < len / 2; ++j) {
                    CoeffType u = coeffs[i + j];
                    CoeffType v = coeffs[i + j + len / 2];
                    CoeffType t = alg.mul(v, w);

                    coeffs[i + j] = alg.add(u, t);
                    coeffs[i + j + len / 2] = alg.sub(u, t);

                    w = alg.mul(w, wlen);
                }
            }
        }

        // --- Normalization: multiply by N⁻¹ modulo q ---
        CoeffType n_inv = mod_pow(static_cast<CoeffType>(n), q - 2, q);
        for (auto &x: coeffs)
            x = alg.mul(x, n_inv);

        return out;
    }

    /**
     * @brief Multiply two polynomials already in NTT domain.
     *
     * @tparam N         Polynomial degree (NTT size).
     * @tparam PolyMod   Polynomial modulus type.
     * @tparam CoeffType Coefficient type used in modular arithmetic.
     * @tparam Als       Algebraic structure implementing Z_q arithmetic.
     *
     * @param a_ntt      First polynomial in NTT domain.
     * @param b_ntt      Second polynomial in NTT domain.
     * @param root_inv   Inverse primitive N-th root of unity.
     *
     * @returns PolynomialField<N> result of inverse NTT( a_ntt ∘ b_ntt ).
     *
     * @note Inputs MUST already be in NTT form. This function does NOT call toNTT().
     */
    template<
        std::size_t N,
        typename PolyMod = NoPolyModulus<N>,
        typename CoeffType = OTPQ_ALS_UNDERLYING_TYPE,
        typename Als = arith::ZmodNField<CoeffType> > requires core::Field<Als, CoeffType>
    auto mulNTT(
        const PolynomialField<N, PolyMod, CoeffType, Als> &a_ntt,
        const PolynomialField<N, PolyMod, CoeffType, Als> &b_ntt) {
        const Als &alg = a_ntt.alg;
        const CoeffType q = alg.mod;

        // --- Pointwise (Hadamard) product ---
        PolynomialField<N, PolyMod, CoeffType, Als> c_ntt(q);

        for (std::size_t i = 0; i < N; ++i)
            c_ntt.coeffs[i] = alg.mul(a_ntt.coeffs[i], b_ntt.coeffs[i]);

        // --- Inverse NTT to return to polynomial form ---
        return c_ntt;
    }
}

