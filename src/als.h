#ifndef OTPQMATH_ALS_H
#define OTPQMATH_ALS_H

#include <concepts>

/**
 * @file als.h
 * @brief Algebraic structure concepts for rings and fields.
 *
 * This header defines compile-time C++23 concepts that formalize
 * algebraic structures used throughout OTPQ-Math:
 *
 *   - otpq::math::als::core::Ring
 *   - otpq::math::als::core::Field
 *
 * These concepts allow generic code to constrain template parameters
 * so they behave like mathematical rings or fields, ensuring that
 * operations such as addition, multiplication, negation, and inversion
 * are available and type-correct.
 *
 * The concepts describe *interfaces* (signatures)–they do not
 * perform algebraic validation (associativity, distributivity, etc.).
 *
 * Intended usage:
 * ----------------
 *   template<Ring<int> R>
 *   void foo(const R& ring);
 *
 *   template<Field<double> F>
 *   void bar(const F& field);
 *
 * A model of Ring or Field is typically a lightweight stateless object
 * (e.g., a Zmod<N> structure) that defines element-wise operations.
 */

namespace otpq::math::als::core {
    // -------------------------------------------------------------------------
    //  Ring Concept
    // -------------------------------------------------------------------------
    /**
     * @brief Concept representing a mathematical ring.
     *
     * A type R models a ring if it provides the following element-wise operations:
     *
     *   - addition:        ring.add(a, b)
     *   - subtraction:     ring.sub(a, b)
     *   - multiplication:  ring.mul(a, b)
     *   - additive inverse: ring.neg(a)
     *   - additive identity: ring.zero()
     *
     * All operations must return values of the element type `Elem`.
     *
     * This concept intentionally stays close to the algebraic definition of a
     * (not necessarily commutative) ring with unity, except that the multiplicative
     * identity is not required. If a multiplicative identity is needed, use Field.
     */
    template<typename R, typename Elem>
    concept Ring =
            requires(const R &ring, Elem a, Elem b)
            {
                /// a + b
                { ring.add(a, b) } -> std::same_as<Elem>;

                /// a - b
                { ring.sub(a, b) } -> std::same_as<Elem>;

                /// a * b
                { ring.mul(a, b) } -> std::same_as<Elem>;

                /// -a   (additive inverse)
                { ring.neg(a) } -> std::same_as<Elem>;

                /// 0    (additive identity)
                { ring.zero() } -> std::same_as<Elem>;
            };


    /**
     * @brief Concept representing a mathematical field.
     *
     * Field extends Ring by additionally requiring:
     *
     *   - multiplicative inverse:  fld.inv(a)
     *   - division:                fld.div(a, b)
     *   - multiplicative identity: fld.one()
     *
     * Every Field is a Ring, but not every Ring is a Field.
     *
     * This maps to the algebraic notion of a commutative unital field
     * (e.g., ℚ, ℝ, finite fields ℤ/pℤ with p prime).
     *
     * NOTE:
     *   Implementers must ensure inv(a) is valid only when a ≠ 0.
     *   Behavior on a = 0 is undefined unless documented otherwise.
     */
    template<typename F, typename Elem>
    concept Field =
            Ring<F, Elem> &&
            requires(const F &fld, Elem a, Elem b)
            {
                /// a⁻¹ (multiplicative inverse)
                { fld.inv(a) } -> std::same_as<Elem>;

                /// a / b
                { fld.div(a, b) } -> std::same_as<Elem>;

                /// 1 (multiplicative identity)
                { fld.one() } -> std::same_as<Elem>;
            };
}


namespace otpq::math::als {
    /**
     * @brief Ring of integers modulo `mod`:  Z / modZ
     *
     * Models the Ring concept for any integral element type `Elem`.
     *
     * Example:
     *     ZmodRing<int> Z7{7};
     *     Z7.add(3, 6);   // 2
     *     Z7.mul(4, 5);   // 6
     */
    template<std::integral Elem>
    struct ZmodRing {
        Elem mod{}; ///< modulus defining the quotient ring

        constexpr explicit ZmodRing(Elem m) noexcept
            : mod{m} {
            // It is allowed to initialize with mod <= 0,
            // but operations with invalid mod are UB.
            // You may enforce (m > 0) with assert or contract if desired.
        }

        // -------------------------------
        // Normalization helper
        // -------------------------------
        [[nodiscard]] constexpr Elem normalize(Elem x) const noexcept {
            x %= mod;
            if (x < 0) x += mod;
            return x;
        }

        // -------------------------------
        // Ring operations
        // -------------------------------
        [[nodiscard]] constexpr Elem zero() const noexcept {
            return Elem{0};
        }

        [[nodiscard]] constexpr Elem neg(Elem a) const noexcept {
            return normalize(-a);
        }

        [[nodiscard]] constexpr Elem add(Elem a, Elem b) const noexcept {
            return normalize(a + b);
        }

        [[nodiscard]] constexpr Elem sub(Elem a, Elem b) const noexcept {
            return normalize(a - b);
        }

        [[nodiscard]] constexpr Elem mul(Elem a, Elem b) const noexcept {
            // safe for any integral: promote to widest type
            using Wide = std::conditional_t<
                (sizeof(Elem) < sizeof(long long)),
                long long,
                __int128_t>;

            Wide r = static_cast<Wide>(a) * static_cast<Wide>(b);
            return normalize(static_cast<Elem>(r % mod));
        }
    };
}


#endif
