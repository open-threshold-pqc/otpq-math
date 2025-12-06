#pragma once

#include <type_traits>
#include <concepts>
#include <stdexcept>

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
 *   - ZmodNField: finite field 𝔽ₙ for prime n
 *   - ComplexField: field of complex numbers
 */

namespace otpq::math::als::ring {
    /**
     * @brief Ring of integers modulo `mod` (ℤ / modℤ).
     *
     * This type models the Ring concept for any integral element type.
     * It supports modular addition, subtraction, multiplication, and
     * negation. Division and multiplicative inverses are *not* provided,
     * since the modulus is not required to be prime.
     *
     * Example:
     * @code
     *   ModRing<int> R{7};   // or simply: ModRing<> R{7}
     *   R.add(3, 6);   // → 2
     *   R.mul(4, 5);   // → 6
     * @endcode
     */
    template<std::integral Elem = OTPQ_ALS_UNDERLYING_TYPE>
    struct ModRing {
        /// @brief The modulus defining ℤ/modℤ.
        Elem mod{};

        /**
         * @brief Constructs the ring with the given modulus.
         * @throws std::invalid_argument if `mod <= 0`.
         */
        explicit constexpr ModRing(Elem m)
            : mod{m} {
            if (mod <= 0)
                throw std::invalid_argument("[ModRing] modulus must be positive");
        }


        /**
         * @brief Returns a unique algebra identifier.
         *
         * Since the modulus fully determines the algebra ℤ / modℤ,
         * this returns the modulus encoded as Elem.
         */
        [[nodiscard]] constexpr Elem id() const noexcept { return mod; }


        /**
         * @brief Reduces an integer into canonical representative `[0, mod)`.
         */
        [[nodiscard]] constexpr Elem normalize(Elem x) const noexcept {
            x %= mod;
            if (x < 0) x += mod;
            return x;
        }

        /// @brief Additive identity (0 mod n).
        [[nodiscard]] constexpr Elem zero() const noexcept { return Elem{0}; }

        /// @brief Additive inverse (-a mod n).
        [[nodiscard]] constexpr Elem neg(Elem a) const noexcept {
            return normalize(-a);
        }

        /// @brief Addition modulo n.
        [[nodiscard]] constexpr Elem add(Elem a, Elem b) const noexcept {
            return normalize(a + b);
        }

        /// @brief Subtraction modulo n.
        [[nodiscard]] constexpr Elem sub(Elem a, Elem b) const noexcept {
            return normalize(a - b);
        }

        /**
         * @brief Multiplication modulo n.
         *
         * Uses a wider intermediate type to avoid overflow when possible.
         */
        [[nodiscard]] constexpr Elem mul(Elem a, Elem b) const noexcept {
            using Wide =
                    std::conditional_t<(sizeof(Elem) < sizeof(long long)),
                        long long,
                        __int128_t>;

            Wide r = static_cast<Wide>(a) * static_cast<Wide>(b);
            return normalize(static_cast<Elem>(r % mod));
        }
    };
}
