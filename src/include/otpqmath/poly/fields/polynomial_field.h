#pragma once

#include <array>
#include <cstddef>
#include <vector>
#include <span>

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
    struct AntiCyclotomicPolyModulus {
        template<typename Coeff, typename Alg>
        [[nodiscard]] std::array<Coeff, N>
        reduce(const std::vector<Coeff> &v, const Alg &alg) const {
            std::array<Coeff, N> out{};
            out.fill(alg.zero());

            for (std::size_t i = 0; i < v.size(); ++i) {
                std::size_t idx = i % N;
                Coeff val = (i < N) ? v[i] : alg.neg(v[i]);
                out[idx] = alg.add(out[idx], val);
            }

            return out;
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
    struct CyclotomicPolyModulus {
        template<typename Coeff, typename Alg>
        [[nodiscard]] std::array<Coeff, N>
        reduce(const std::vector<Coeff> &v, const Alg &alg) const {
            std::array<Coeff, N> out{};
            out.fill(alg.zero());

            for (std::size_t i = 0; i < v.size(); ++i) {
                out[i % N] = alg.add(out[i % N], v[i]);
            }

            return out;
        }

        template<typename Coeff, typename Alg>
        constexpr void reduce(const std::array<Coeff, N> &, const Alg &) const {
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
        template<typename Coeff, typename Alg>
        [[nodiscard]] std::array<Coeff, N>
        reduce(const std::vector<Coeff> &v, const Alg &alg) const {
            std::array<Coeff, N> out{};
            out.fill(alg.zero());

            //todo if the mult does not fit
            if (v.size() > N)
                throw std::runtime_error(
                    std::format("[PolynomialField] * result with no poly modulus does not fit {} coefficients", N));

            for (std::size_t i = 0; i < std::min(N, v.size()); ++i)
                out[i] = v[i];

            return out;
        }

        template<typename Coeff, typename Alg>
        constexpr void reduce(const std::array<Coeff, N> &, const Alg &) const {
        }
    };

    /**
     * @class PolynomialField
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
     * @tparam PolyMod   Modulus policy implementing reduction
     * @tparam CoeffType Coefficient type used in storage
     * @tparam Als       Algebra defining coefficient operations
     *
     * Requirements:
     *   • Als must satisfy the `core::Field` concept.
     *   • All constructors explicitly require an `Als` instance.
     */
    template<
        std::size_t N,
        typename PolyMod = NoPolyModulus<N>,
        typename CoeffType = OTPQ_ALS_UNDERLYING_TYPE,
        typename Als = field::PrimeField<CoeffType> >
        requires core::Field<Als, CoeffType>
    class PolynomialField : public Polynomial<N, CoeffType> {
    public:
        using Base = Polynomial<N, CoeffType>;
        using Base::coeffs_;

        /**
         * @brief Default constructor is disabled.
         *
         * An algebra instance must always be explicitly provided so that
         * coefficient arithmetic is well-defined.
         */
        PolynomialField() = delete;

        /**
         * @brief Constructs the zero polynomial with explicit algebra and modulus.
         *
         * @param a  Coefficient algebra.
         * @param m  Polynomial modulus policy.
         */
        constexpr explicit PolynomialField(const Als &a,
                                           const PolyMod &m = PolyMod())
            : Base(), als_(a), polyModulus_(m) {
        }

        /**
         * @brief Constructs from an array of N coefficients and reduces.
         *
         * @param arr Coefficient array.
         * @param a   Coefficient algebra.
         * @param m   Polynomial modulus policy.
         */
        constexpr explicit PolynomialField(const std::array<CoeffType, N> &arr,
                                           const Als &a,
                                           const PolyMod &m = PolyMod())
            : Base(arr), als_(a), polyModulus_(m) {
            reduce();
        }

        /**
         * @brief Constructs from a span of exactly N coefficients and reduces.
         *
         * @param s   Span of N coefficients.
         * @param a   Coefficient algebra.
         * @param m   Polynomial modulus policy.
         */
        constexpr explicit PolynomialField(std::span<const CoeffType, N> s,
                                           const Als &a,
                                           const PolyMod &m = PolyMod())
            : Base(s), als_(a), polyModulus_(m) {
            reduce();
        }

        /**
         * @brief Constructs from an initializer list (truncated or padded) and reduces.
         *
         * @param list Initializer list of coefficients.
         * @param a    Coefficient algebra.
         * @param m    Polynomial modulus policy.
         */
        constexpr explicit PolynomialField(std::initializer_list<CoeffType> list,
                                           const Als &a,
                                           const PolyMod &m = PolyMod())
            : Base(list), als_(a), polyModulus_(m) {
            reduce();
        }

        /**
         * @brief Constructs constant polynomial p(x) = constant_value and reduces.
         *
         * @param constant_value Constant term.
         * @param a              Coefficient algebra.
         * @param m              Polynomial modulus policy.
         */
        constexpr explicit PolynomialField(CoeffType constant_value,
                                           const Als &a,
                                           const PolyMod &m = PolyMod())
            : Base(constant_value), als_(a), polyModulus_(m) {
            reduce();
        }

        /**
         * @brief Construct from an arbitrary container (truncated/padded).
         *
         * Copies up to N coefficients from the container; excess elements are
         * ignored and missing elements are zero. Supports any iterable container
         * except another PolynomialField.
         */
        template<typename Container>
            requires (!std::is_same_v<std::decay_t<Container>, PolynomialField>)
        constexpr explicit PolynomialField(const Container &c,
                                           const Als &a,
                                           const PolyMod &m = PolyMod())
            : Base(), als_(a), polyModulus_(m) {
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
        [[nodiscard]] constexpr const Als &algebra() const noexcept {
            return als_;
        }

        /**
         * @brief Returns the polynomial modulus reduction policy.
         */
        [[nodiscard]] constexpr const PolyMod &modulus() const noexcept {
            return polyModulus_;
        }

        /**
         * @brief Returns true if the polynomial is identically zero.
         *
         * Zero is defined as having degree 0 and its constant term equal to zero
         * in the coefficient algebra.
         */
        [[nodiscard]] constexpr bool is_zero() const {
            return this->degree() == 0 && als_.is_zero(coeffs_[0]);
        }

        /**
         * @brief Addition in Alg[x] / M(x).
         *
         * @param a Left operand.
         * @param b Right operand.
         * @return  Reduced sum.
         */
        friend constexpr PolynomialField
        operator+(const PolynomialField &a, const PolynomialField &b) {
            PolynomialField r(a.algebra(), a.modulus());
            for (std::size_t i = 0; i < N; ++i)
                r.coeffs_[i] = a.als_.add(a.coeffs_[i], b.coeffs_[i]);
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
        friend constexpr PolynomialField
        operator-(const PolynomialField &a, const PolynomialField &b) {
            PolynomialField r(a.als_, a.polyModulus_);
            for (std::size_t i = 0; i < N; ++i)
                r.coeffs_[i] = a.als_.sub(a.coeffs_[i], b.coeffs_[i]);
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
        friend PolynomialField operator*(const PolynomialField &a,
                                         const PolynomialField &b) {
            const std::size_t da = a.degree();
            const std::size_t db = b.degree();

            std::vector<CoeffType> tmp(da + db + 1, a.als_.zero());

            for (std::size_t i = 0; i <= da; ++i)
                for (std::size_t j = 0; j <= db; ++j)
                    tmp[i + j] = a.als_.add(tmp[i + j],
                                            a.als_.mul(a.coeffs_[i], b.coeffs_[j]));

            auto reduced = a.modulus().reduce(tmp, a.als_);
            return PolynomialField(reduced, a.als_, a.polyModulus_);
        }

    private:
        Als als_; /**< Coefficient algebra (field or ring). */
        PolyMod polyModulus_; /**< Polynomial modulus reduction policy. */

        /**
         * @brief Applies polynomial modulus reduction to the coefficient array.
         */
        constexpr void reduce() { polyModulus_.reduce(coeffs_, als_); }
    };

    /**
     * @brief Cyclotomic polynomial ring Alg[x] / (x^N - 1).
     */
    template<std::size_t N,
        typename CoeffType = OTPQ_ALS_UNDERLYING_TYPE,
        typename Als = field::PrimeField<CoeffType> >
    using PolynomialFieldCyclotomic =
    PolynomialField<N, CyclotomicPolyModulus<N>, CoeffType, Als>;
}
