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
     *   - normalize(x)   → canonical coefficient reduction (e.g., modular normalization)
     *
     * Additionally, every algebraic ring must provide a unique identifier:
     *
     *   - id()           → a runtime identifier describing the algebra instance
     *
     * The identifier **uniquely distinguishes algebraic structures**, and is used
     * to ensure compatibility when performing polynomial-ring operations. For
     * modular rings such as ℤ_q, id() will typically encode or hash the modulus q.
     * For parameter-free rings (e.g., ℝ, ℂ), id() may return a fixed constant
     * representing that algebra.
     *
     * Two ring instances R₁ and R₂ are considered compatible iff:
     *
     *        R₁.id() == R₂.id()
     *
     * This concept specifies the interface only; it does not enforce
     * algebraic laws such as associativity or distributivity.
     */
    template<typename R, typename Elem>
    concept Ring =
            std::equality_comparable<Elem> &&
            requires(const R &ring, Elem a, Elem b)
            {
                { ring.add(a, b) } -> std::same_as<Elem>;
                { ring.sub(a, b) } -> std::same_as<Elem>;
                { ring.mul(a, b) } -> std::same_as<Elem>;
                { ring.neg(a) } -> std::same_as<Elem>;
                { ring.zero() } -> std::same_as<Elem>;
                { ring.normalize(a) } -> std::same_as<Elem>; // Maybe a no-op
                { ring.id() } -> std::same_as<Elem>;
            };

    /**
     * @brief Concept representing a mathematical Field over type Elem.
     *
     * A Field extends Ring by requiring:
     *
     *   one()    → multiplicative identity
     *   inv(a)   → multiplicative inverse (a ≠ 0)
     *   div(a,b) → a / b
     *
     * This concept is satisfied by:
     *   - ℤ/pℤ for prime p (ZmodNField)
     *   - ℝ, ℚ (if someone implements)
     *   - ℂ (ComplexAls)
     *
     * Any valid ALS for a Field must document the behavior of inv(a)
     * when a == 0 (undefined in general).
     */
    template<typename F, typename Elem>
    concept Field = Ring<F, Elem> && requires(const F &fld, Elem a, Elem b)
    {
        { fld.inv(a) } -> std::same_as<Elem>;
        { fld.div(a, b) } -> std::same_as<Elem>;
        { fld.one() } -> std::same_as<Elem>;
    };


    /**
     * @brief Concept: ANY algebraic structure is either Field or Ring
    *
    */
    template<typename Als, typename Elem>
    concept AlgebraicStructure =
            Ring<Als, Elem> || Field<Als, Elem>;
}

// //todo to be thought
// if constexpr (Als::supports_simd) {
//     simd_add(r.coeffs_, a.coeffs_, b.coeffs_, als_.q());
// } else {
//     // fallback
//     for (i...)
//         r[i] = als_.add(a[i], b[i]);
// }
