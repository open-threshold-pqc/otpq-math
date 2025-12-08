#pragma once

#include <type_traits>
#include <concepts>
#include <stdexcept>
#include <otpqmath/utils/math_utils.h>

/**
 * @brief Integer rings and finite fields modulo `n`.
 *
 * This namespace provides modular arithmetic structures based on
 * integer types. These models satisfy the algebraic concepts defined
 * in `otpq::math::als::core` such as Ring and Field, and are used
 * throughout OTPQ-Math for polynomial arithmetic and number theory.
 *
 * Included structures:
 *   - ModRing:  additive/multiplicative ring ℤ/nℤ
 *   - PrimeField: finite field 𝔽ₙ for prime n
 *   - ComplexField: field of complex numbers
 */

namespace otpq::math::als::field::number {
    /**
     * @brief Finite field of integers modulo prime `mod` (𝔽_mod).
     *
     * This type models the Field concept when `mod` is prime. It supports
     * modular addition, subtraction, multiplication, negation, division,
     * and multiplicative inversion.
     *
     * Example:
     * @code
     *   PrimeField<int> F7{7};   // or simply:   PrimeField<> F7{7}
     *   F7.add(3, 6);   // → 2
     *   F7.mul(4, 5);   // → 6
     *   F7.div(3, 4);   // → 5  (since 4^{-1} = 2 mod 7, and 3*2 = 6 ≡ 6)
     * @endcode
     */
    template<std::integral Elem = OTPQ_ALS_UNDERLYING_TYPE>
    struct PrimeField {
        /// @brief Prime modulus defining the finite field 𝔽ₘₒd.
        Elem mod{};

        /**
         * @brief Constructs a finite field with prime modulus.
         * @throws std::invalid_argument if modulus ≤ 1 or not prime.
         */
        explicit constexpr PrimeField(Elem m)
            : mod{m} {
            if (mod <= 1)
                throw std::invalid_argument("[PrimeField] modulus must be > 1");

            if (!utils::math::is_prime(mod))
                throw std::invalid_argument("[PrimeField] modulus must be prime");
        }

        /**
         * @brief Returns a unique algebra identifier.
         *
         * Since the modulus fully determines the algebra ℤ / modℤ,
         * this returns the modulus encoded as Elem.
         */
        [[nodiscard]] constexpr Elem id() const noexcept { return mod; }


        /// @brief Reduces an integer into canonical representative `[0, mod)`.
        [[nodiscard]] constexpr Elem normalize(Elem x) const noexcept {
            x %= mod;
            if (x < 0) x += mod;
            return x;
        }

        /// @brief Additive identity.
        [[nodiscard]] constexpr Elem zero() const noexcept { return Elem{0}; }

        /// @brief Multiplicative identity.
        [[nodiscard]] constexpr Elem one() const noexcept { return Elem{1}; }

        /// @brief Additive inverse (-a mod p).
        [[nodiscard]] constexpr Elem neg(Elem a) const noexcept {
            return normalize(-a);
        }

        /// @brief Addition modulo p.
        [[nodiscard]] constexpr Elem add(Elem a, Elem b) const noexcept {
            return normalize(a + b);
        }

        /// @brief Subtraction modulo p.
        [[nodiscard]] constexpr Elem sub(Elem a, Elem b) const noexcept {
            return normalize(a - b);
        }

        /**
         * @brief Multiplication modulo p.
         *
         * Promotes to a wide intermediate type to avoid overflow.
         */
        [[nodiscard]] constexpr Elem mul(Elem a, Elem b) const noexcept {
            using Wide =
                    std::conditional_t<(sizeof(Elem) < sizeof(long long)),
                        long long,
                        __int128_t>;

            Wide wa = static_cast<Wide>(a);
            Wide wb = static_cast<Wide>(b);

            Wide r = (wa * wb) % static_cast<Wide>(mod);
            return normalize(static_cast<Elem>(r));
        }

        /**
         * @brief Multiplicative inverse via extended Euclidean algorithm.
         * @throws std::runtime_error if the element is not invertible.
         */
        [[nodiscard]] constexpr Elem inv(Elem a) const {
            Elem t = 0, newt = 1;
            Elem r = mod, newr = normalize(a);

            while (newr != 0) {
                Elem q = r / newr;

                Elem old_t = t;
                t = newt;
                newt = old_t - q * newt;

                Elem old_r = r;
                r = newr;
                newr = old_r - q * newr;
            }

            if (r > 1)
                throw std::runtime_error("PrimeField: element not invertible");

            if (t < 0)
                t += mod;

            return normalize(t);
        }

        /// @brief Division modulo p: `a / b = a * inv(b) mod p`.
        [[nodiscard]] constexpr Elem div(Elem a, Elem b) const {
            return mul(a, inv(b));
        }
    };
}
