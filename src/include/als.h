#pragma once

#include <concepts>


/**
 * @brief Concepts describing algebraic structures (als) interfaces such as Ring and Field.
 *
 * Here, we define pure compile-time concepts specifying
 * the required operations for algebraic structures. These concepts are used
 * throughout the library to ensure that template parameters conform to the
 * expected mathematical interfaces.
 *
 * No concrete implementations are provided here, only the abstract
 * requirements for:
 *   - Ring (add, sub, mul, neg, zero)
 *   - Field (Ring + inv, div, one)
 *
 * We note that other properties such as associativity, distribution under multiplication,
 * ring commutativity (if exists), etc. have not been considered for ease of implementation.
 *
 */
namespace otpq::math::als::core {
    /**
     * @brief Concept representing an algebraic ring over element type Elem.
     *
     * A type R models a ring if it provides element-wise operations
     * corresponding to the algebraic structure of a (not necessarily
     * commutative) ring with additive identity:
     *
     *   - add(a, b) → a + b
     *   - sub(a, b) → a - b
     *   - mul(a, b) → a * b
     *   - neg(a)    → -a
     *   - zero()    → additive identity
     *
     * All operations must return Elem.
     *
     * This concept specifies the interface only; it does not enforce
     * algebraic laws such as associativity or distributivity.
     */
    template<typename R, typename Elem>
    concept Ring = requires(const R &ring, Elem a, Elem b)
    {
        { ring.add(a, b) } -> std::same_as<Elem>;
        { ring.sub(a, b) } -> std::same_as<Elem>;
        { ring.mul(a, b) } -> std::same_as<Elem>;
        { ring.neg(a) } -> std::same_as<Elem>;
        { ring.zero() } -> std::same_as<Elem>;
    };


    /**
     * @brief Concept representing a mathematical field over element type Elem.
     *
     * A Field extends Ring by requiring multiplicative identity, division,
     * and inversion:
     *
     *   - one()    → multiplicative identity
     *   - inv(a)   → multiplicative inverse (a ≠ 0)
     *   - div(a,b) → a / b
     *
     * Fields model the algebraic notion of a commutative field such as ℚ,
     * ℝ, or finite fields ℤ/pℤ (with p prime).
     *
     * Implementations must document the behavior of inv(a) when a = 0,
     * as this is undefined in general.
     */
    template<typename F, typename Elem>
    concept Field = Ring<F, Elem> && requires(const F &fld, Elem a, Elem b)
    {
        { fld.inv(a) } -> std::same_as<Elem>;
        { fld.div(a, b) } -> std::same_as<Elem>;
        { fld.one() } -> std::same_as<Elem>;
    };
}
