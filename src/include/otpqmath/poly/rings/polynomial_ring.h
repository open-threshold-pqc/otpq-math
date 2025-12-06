#pragma once

#include <array>
#include <cstddef>
#include <vector>
#include <span>
#include <algorithm>

#include <otpqmath/algebra/als.h>
#include <otpqmath/poly/polynomial.h>

namespace otpq::math::als::polynomial {
    /**
     * @brief Reduction modulo (xⁿ + 1).
     *
     * Implements the relation
     *
     *      xⁿ = -1   →   xⁿ⁺ᵏ = -xᵏ
     *
     * All coefficients beyond index N−1 are folded back with a negation.
     * Commonly used in lattice cryptography (e.g., NTRU, RLWE variants).
     */
    template<std::size_t N>
    struct PolyAntiCyclotomicModulus {
        /**
         * @brief Reduces a variable-length coefficient vector modulo xⁿ + 1.
         *
         * Used for operations (e.g., multiplication) that produce more than N
         * coefficients. Performs structural negacyclic folding and returns an
         * N-coefficient canonical representative.
         */
        template<typename CoeffType, typename CoeffAls>
            requires core::AlgebraicStructure<CoeffAls, CoeffType>
        [[nodiscard]] std::array<CoeffType, N>
        reduce(const std::vector<CoeffType> &coeffs, const CoeffAls &coeffAls) const noexcept {
            std::array<CoeffType, N> reducedCoeffs{};
            reducedCoeffs.fill(coeffAls.zero());

            for (std::size_t i = 0; i < coeffs.size(); ++i) {
                std::size_t idx = i % N;
                bool negate = ((i / N) & 1) != 0;
                CoeffType val = negate ? coeffAls.neg(coeffs[i]) : coeffs[i];
                reducedCoeffs[idx] = coeffAls.add(reducedCoeffs[idx], val);
            }

            return reducedCoeffs;
        }

        /**
         * @brief Normalizes an N-coefficient polynomial in-place.
         *
         * Used for operations that preserve size N (e.g., initialization, addition,
         * subtraction). Applies coefficient-level normalization only.
         */
        template<typename CoeffType, typename CoeffAls>
            requires core::AlgebraicStructure<CoeffAls, CoeffType>
        constexpr void reduce(std::array<CoeffType, N> &coeffs, const CoeffAls &coeffAls) const noexcept {
            for (auto &coeff: coeffs)
                coeff = coeffAls.normalize(coeff);
        }
    };

    /**
     * @brief Reduction modulo (xⁿ − 1).
     *
     * Wrap-around relation:
     *
     *      xⁿ = 1   →   xⁿ⁺ᵏ = xᵏ
     *
     * All coefficients are accumulated by index modulo N.
     */
    template<std::size_t N>
    struct PolyCyclotomicModulus {
        /**
         * @brief Reduces a variable-length coefficient vector modulo xⁿ - 1.
         *
         * Used for operations (e.g., multiplication) that produce more than N
         * coefficients. Performs structural cyclic folding and returns an
         * N-coefficient canonical representative.
         */
        template<typename CoeffType, typename CoeffAls>
            requires core::AlgebraicStructure<CoeffAls, CoeffType>
        [[nodiscard]] std::array<CoeffType, N>
        reduce(const std::vector<CoeffType> &coeffs, const CoeffAls &coeffAls) const noexcept {
            std::array<CoeffType, N> reducedCoeffs{};
            reducedCoeffs.fill(coeffAls.zero());

            for (std::size_t i = 0; i < coeffs.size(); ++i)
                reducedCoeffs[i % N] = coeffAls.add(reducedCoeffs[i % N], coeffs[i]);

            return reducedCoeffs;
        }

        /**
         * @brief Normalizes an N-coefficient polynomial in-place.
         *
         * Used for operations that preserve size N (e.g., initialization, addition,
         * subtraction). Applies coefficient-level normalization only.
         */
        template<typename CoeffType, typename CoeffAls>
            requires core::AlgebraicStructure<CoeffAls, CoeffType>
        constexpr void reduce(std::array<CoeffType, N> &coeffs, const CoeffAls &coeffAls) const noexcept {
            for (auto &coeff: coeffs)
                coeff = coeffAls.normalize(coeff);
        }
    };

    /**
     * @brief No modulus (raw polynomial truncation).
     *
     * Copies up to N coefficients from the input vector.
     * Useful for representing polynomials without a quotient structure.
     */
    template<std::size_t N>
    struct NoPolyModulus {
        template<typename CoeffType, typename Alg>
        [[nodiscard]] std::array<CoeffType, N>
        reduce(const std::vector<CoeffType> &coeffs, const Alg &alg) const {
            std::array<CoeffType, N> reducedCoeffs{};
            reducedCoeffs.fill(alg.zero());

            //todo if the mult does not fit
            if (coeffs.size() > N)
                throw std::runtime_error(
                    std::format("[PolynomialRing] * result with no poly modulus does not fit {} coefficients", N));

            for (std::size_t i = 0; i < std::min(N, coeffs.size()); ++i)
                reducedCoeffs[i] = coeffs[i];

            return reducedCoeffs;
        }

        /**
         * @brief Normalizes an N-coefficient polynomial in-place.
         *
         * Used for operations that preserve size N (e.g., initialization, addition,
         * subtraction). Applies coefficient-level normalization only.
         */
        template<typename CoeffType, typename CoeffAls>
            requires core::AlgebraicStructure<CoeffAls, CoeffType>
        constexpr void reduce(std::array<CoeffType, N> &coeffs, const CoeffAls &coeffAls) const noexcept {
            for (auto &coeff: coeffs)
                coeff = coeffAls.normalize(coeff);
        }
    };
}


namespace otpq::math::als::polynomial {
    /**
     * @class PolynomialRing
     * @brief Polynomial element of the quotient ring
     *
     *        Alg[x] / M(x)
     *
     * where:
     *   • Alg     – an algebraic structure providing coefficient arithmetic
     *               (e.g., ZmodNField)
     *   • M(x)    – a polynomial modulus determined by `PolyMod`
     *   • N       – static coefficient capacity (degree < N)
     *
     * This class extends `Polynomial<N>` with algebraic operations:
     * addition, subtraction, multiplication, and reduction modulo M(x).
     *
     * Coefficients are stored in the base class `Polynomial<N>`.  The `PolyMod`
     * policy provides polynomial reduction (e.g., xⁿ − 1, xⁿ + 1, no modulus).
     *
     * @tparam N         Polynomial capacity (static)
     * @tparam PolyModulusPolicy   Modulus policy implementing reduction
     * @tparam CoeffType Coefficient type used in storage
     * @tparam CoeffAls       Algebra defining coefficient operations
     *
     * Requirements:
     *   • Als must satisfy the `core::Field` concept.
     *   • All constructors explicitly require an `Als` instance.
     */
    template<
        std::size_t N,
        typename PolyModulusPolicy = NoPolyModulus<N>,
        typename CoeffType = OTPQ_ALS_UNDERLYING_TYPE,
        typename CoeffAls = field::PrimeField<CoeffType> >
        requires core::AlgebraicStructure<CoeffAls, CoeffType>
    class PolynomialRing : public Polynomial<N, CoeffType> {
    public:
        using Base = Polynomial<N, CoeffType>;
        using Base::coeffs_;

        /**
         * @brief Default constructor is disabled.
         *
         * An algebra instance must always be explicitly provided so that
         * coefficient arithmetic is well-defined.
         */
        PolynomialRing() = delete;

        /**
         * @brief Constructs the zero polynomial with explicit algebra and modulus.
         *
         * @param a  Coefficient algebra.
         * @param m  Polynomial modulus policy.
         */
        constexpr explicit PolynomialRing(const CoeffAls &a,
                                          const PolyModulusPolicy &m = PolyModulusPolicy())
            : Base(), coeffAls_(a), polyModPolicy_(m) {
        }

        /**
         * @brief Constructs from an array of N coefficients and reduces.
         *
         * @param arr Coefficient array.
         * @param a   Coefficient algebra.
         * @param m   Polynomial modulus policy.
         */
        constexpr explicit PolynomialRing(const std::array<CoeffType, N> &arr,
                                          const CoeffAls &a,
                                          const PolyModulusPolicy &m = PolyModulusPolicy())
            : Base(arr), coeffAls_(a), polyModPolicy_(m) {
            reduce();
        }

        /**
         * @brief Constructs from a span of exactly N coefficients and reduces.
         *
         * @param s   Span of N coefficients.
         * @param a   Coefficient algebra.
         * @param m   Polynomial modulus policy.
         */
        constexpr explicit PolynomialRing(std::span<const CoeffType, N> s,
                                          const CoeffAls &a,
                                          const PolyModulusPolicy &m = PolyModulusPolicy())
            : Base(s), coeffAls_(a), polyModPolicy_(m) {
            reduce();
        }

        /**
         * @brief Constructs from an initializer list (truncated or padded) and reduces.
         *
         * @param list Initializer list of coefficients.
         * @param a    Coefficient algebra.
         * @param m    Polynomial modulus policy.
         */
        constexpr explicit PolynomialRing(std::initializer_list<CoeffType> list,
                                          const CoeffAls &a,
                                          const PolyModulusPolicy &m = PolyModulusPolicy())
            : Base(list), coeffAls_(a), polyModPolicy_(m) {
            reduce();
        }

        /**
         * @brief Constructs constant polynomial p(x) = constant_value and reduces.
         *
         * @param constant_value Constant term.
         * @param a              Coefficient algebra.
         * @param m              Polynomial modulus policy.
         */
        constexpr explicit PolynomialRing(CoeffType constant_value,
                                          const CoeffAls &a,
                                          const PolyModulusPolicy &m = PolyModulusPolicy())
            : Base(constant_value), coeffAls_(a), polyModPolicy_(m) {
            reduce();
        }

        /**
         * @brief Construct from an arbitrary container (truncated/padded).
         *
         * Copies up to N coefficients from the container; excess elements are
         * ignored and missing elements are zero. Supports any iterable container
         * except another PolynomialRing.
         */
        template<typename Container>
            requires (!std::is_same_v<std::decay_t<Container>, PolynomialRing>)
        constexpr explicit PolynomialRing(const Container &c,
                                          const CoeffAls &a,
                                          const PolyModulusPolicy &m = PolyModulusPolicy())
            : Base(), coeffAls_(a), polyModPolicy_(m) {
            std::size_t i = 0;

            for (auto &&v: c) {
                if (i >= N) break;
                this->coeffs_[i++] = static_cast<CoeffType>(v);
            }
            reduce();
        }

        /**
         * @brief Returns the coefficient algebra used by this polynomial.
         */
        [[nodiscard]] constexpr const CoeffAls &coeff_als() const noexcept {
            return coeffAls_;
        }

        /**
         * @brief Returns the polynomial modulus reduction policy.
         */
        [[nodiscard]] constexpr const PolyModulusPolicy &poly_mod_policy() const noexcept {
            return polyModPolicy_;
        }

        /**
         * @brief Returns true if the polynomial is identically zero.
         *
         * Zero is defined as having degree 0 and its constant term equal to zero
         * in the coefficient algebra.
         */
        [[nodiscard]] constexpr bool is_zero() const {
            return this->degree() == 0 && (coeffs_[0] == coeffAls_.zero());
        }

        /**
         * @brief Addition in Alg[x] / M(x).
         *
         * @param a Left operand.
         * @param b Right operand.
         * @return  Reduced sum.
         */
        friend constexpr PolynomialRing
        operator+(const PolynomialRing &a, const PolynomialRing &b) {
            if (a.coeffAls_.id() != b.coeffAls_.id())
                throw std::runtime_error(
                    "[PolynomialRing] operator+: incompatible coefficient algebra (CoeffAls.id mismatch)");
            PolynomialRing r(a.coeffAls_, a.polyModPolicy_);
            for (std::size_t i = 0; i < N; ++i)
                r.coeffs_[i] = a.coeffAls_.add(a.coeffs_[i], b.coeffs_[i]);
            r.reduce();
            return r;
        }

        /**
         * @brief Subtraction in Alg[x] / M(x).
         *
         * @param a Left operand.
         * @param b Right operand.
         * @return  Reduced difference.
         */
        friend constexpr PolynomialRing
        operator-(const PolynomialRing &a, const PolynomialRing &b) {
            if (a.coeffAls_.id() != b.coeffAls_.id())
                throw std::runtime_error(
                    "[PolynomialRing] operator-: incompatible coefficient algebra (CoeffAls.id mismatch)");

            PolynomialRing r(a.coeffAls_, a.polyModPolicy_);
            for (std::size_t i = 0; i < N; ++i)
                r.coeffs_[i] = a.coeffAls_.sub(a.coeffs_[i], b.coeffs_[i]);
            r.reduce();
            return r;
        }

        /**
         * @brief Multiplication in Alg[x] / M(x) via naive convolution.
         *
         * Performs a full convolution of degrees ≤ 2N−2 and applies the modulus
         * reduction defined by `PolyMod`. The result is a canonical representative
         * modulo M(x).
         *
         * @param a Left operand.
         * @param b Right operand.
         * @return  Reduced product.
         */
        friend PolynomialRing operator*(const PolynomialRing &a,
                                        const PolynomialRing &b) {
            if (a.coeffAls_.id() != b.coeffAls_.id())
                throw std::runtime_error(
                    "[PolynomialRing] operator*: incompatible coefficient algebra (CoeffAls.id mismatch)");

            const std::size_t da = a.degree();
            const std::size_t db = b.degree();

            std::vector<CoeffType> tmp(da + db + 1, a.coeffAls_.zero());

            for (std::size_t i = 0; i <= da; ++i)
                for (std::size_t j = 0; j <= db; ++j)
                    tmp[i + j] = a.coeffAls_.add(tmp[i + j],
                                                 a.coeffAls_.mul(a.coeffs_[i], b.coeffs_[j]));

            auto reduced = a.polyModPolicy_.reduce(tmp, a.coeffAls_);
            return PolynomialRing(reduced, a.coeffAls_, a.polyModPolicy_);
        }

    private:
        CoeffAls coeffAls_; /**< Coefficient algebra (field or ring). */
        PolyModulusPolicy polyModPolicy_; /**< Polynomial modulus reduction policy. */

        /**
         * @brief Applies polynomial modulus reduction to the coefficient array.
         */
        constexpr void reduce() { polyModPolicy_.reduce(coeffs_, coeffAls_); }
    };

    /**
     * @brief Cyclotomic polynomial ring Als[x] / (x^N - 1).
     */
    template<std::size_t N,
        typename CoeffType = OTPQ_ALS_UNDERLYING_TYPE,
        typename Als = field::PrimeField<CoeffType> >
    using PolynomialCyclotomicRing =
    PolynomialRing<N, PolyCyclotomicModulus<N>, CoeffType, Als>;
}
