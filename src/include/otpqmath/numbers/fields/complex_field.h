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
 *   - ZmodNRing:  additive/multiplicative ring ℤ/nℤ
 *   - ZmodNField: finite field 𝔽ₙ for prime n
 *   - ComplexField: field of complex numbers
 */

namespace otpq::math::als::field {
    /**
     * @brief Algebraic structure for complex-number arithmetic (ℂ).
     *
     * This ALS implements the `Field` concept over the element type:
     *
     *      Elem = std::complex<Real>
     *
     * It provides:
     *   - addition, subtraction, negation
     *   - multiplication, division
     *   - additive identity (0)
     *   - multiplicative identity (1)
     *   - multiplicative inverse (1/a)
     *
     * This field is essential for algorithms operating in the Fourier domain,
     * such as:
     *   - FFT-based polynomial multiplication
     *   - Falcon’s recursive FFT splitting (complex 1024-point FFT)
     *   - GPV sampling over ℤⁿ using Gram–Schmidt represented in ℂ
     *
     * Unlike modular arithmetic, ℂ is a true mathematical field:
     *
     *      ∀ a ≠ 0,  ∃ a⁻¹  such that  a · a⁻¹ = 1
     *
     * Therefore it satisfies your `core::Field` concept without modification.
     *
     * @tparam Real  The underlying real type for complex components.
     *               Typically `double` (for FFT precision).
     */
    template<typename Real = double>
    struct ComplexField {
        /// @brief Complex field element type.
        using Elem = std::complex<Real>;

        /**
         * @brief Returns a unique algebra identifier.
         *
         * Since the modulus fully determines the algebra of complex numbers.
         */
        [[nodiscard]] constexpr Elem id() const noexcept { return Elem{0}; }

        /**
         * @brief Returns additive identity 0 + 0i.
         *
         * @return Elem(0,0)
         */
        [[nodiscard]] constexpr Elem zero() const noexcept {
            return Elem{Real{0}, Real{0}};
        }

        /**
         * @brief Returns multiplicative identity 1 + 0i.
         *
         * @return Elem(1,0)
         */
        [[nodiscard]] constexpr Elem one() const noexcept {
            return Elem{Real{1}, Real{0}};
        }

        // -------------------------------------------------------------------------
        // Basic operations (satisfying Ring interface)
        // -------------------------------------------------------------------------

        /**
         * @brief Complex addition: a + b
         */
        [[nodiscard]] constexpr Elem add(const Elem &a, const Elem &b) const noexcept {
            return a + b;
        }

        /**
         * @brief Complex subtraction: a - b
         */
        [[nodiscard]] constexpr Elem sub(const Elem &a, const Elem &b) const noexcept {
            return a - b;
        }

        /**
         * @brief Complex multiplication: a * b
         *
         * Uses built-in std::complex operator*.
         */
        [[nodiscard]] constexpr Elem mul(const Elem &a, const Elem &b) const noexcept {
            return a * b;
        }

        /**
         * @brief Additive inverse: -a
         */
        [[nodiscard]] constexpr Elem neg(const Elem &a) const noexcept {
            return -a;
        }

        // -------------------------------------------------------------------------
        // Field operations (inverse and division)
        // -------------------------------------------------------------------------

        /**
         * @brief Multiplicative inverse: 1 / a
         *
         * @note Undefined if a == 0.
         *
         * Behavior mirrors std::complex<double> semantics. Caller must ensure a ≠ 0.
         */
        [[nodiscard]] constexpr Elem inv(const Elem &a) const noexcept {
            return Elem{Real{1}, Real{0}} / a;
        }

        /**
         * @brief Division: a / b
         *
         * Equivalent to a * inv(b).
         */
        [[nodiscard]] constexpr Elem div(const Elem &a, const Elem &b) const noexcept {
            return a / b;
        }
    };
}
