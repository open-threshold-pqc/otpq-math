#pragma once

#include <array>
#include <cstddef>
#include <vector>

#include <zmod_arithmetic.h>
#include <polynomial.h>
#include <als.h>

namespace otpq::math::als::polynomial {
    /**
     * @brief Reduction mod (x^N + 1) for polynomials of dynamic length.
     *
     * Uses the rule:
     *      x^(N+k) = -x^k
     *
     * Useful in lattice cryptography (e.g., NTT variants, RLWE).
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
                Coeff val = v[i];

                if (i >= N)
                    val = alg.neg(val);

                out[idx] = alg.add(out[idx], val);
            }

            return out;
        }
    };

    /**
     * @brief Reduction mod (x^N - 1).
     *
     * Wrap-around follows:
     *      c[i mod N] += c[i]
     */
    template<std::size_t N>
    struct CyclotomicPolyModulus {
        template<typename Coeff, typename Alg>
        [[nodiscard]] std::array<Coeff, N>
        reduce(const std::vector<Coeff> &v, const Alg &alg) const {
            std::array<Coeff, N> out{};
            out.fill(alg.zero());

            for (std::size_t i = 0; i < v.size(); ++i) {
                std::size_t idx = i % N;
                out[idx] = alg.add(out[idx], v[i]);
            }

            return out;
        }

        template<typename Coeff, typename Alg>
        constexpr void reduce(const std::array<Coeff, N> &, const Alg &) const {
        }
    };

    /**
     * @brief No polynomial modulus — copies up to N coefficients.
     *
     * Used for raw polynomials without quotient structure.
     */
    template<std::size_t N>
    struct NoPolyModulus {
        template<typename Coeff, typename Alg>
        [[nodiscard]] std::array<Coeff, N>
        reduce(const std::vector<Coeff> &v, const Alg &alg) const {
            std::array<Coeff, N> out{};
            out.fill(alg.zero());

            // todo check for fitness !

            for (std::size_t i = 0; i < std::min(N, v.size()); ++i)
                out[i] = v[i];

            return out;
        }

        template<typename Coeff, typename Alg>
        constexpr void reduce(const std::array<Coeff, N> &, const Alg &) const {
        }
    };

    /**
     * @brief Polynomial arithmetic over an algebraic field or ring.
     *
     * Represents polynomials in:
     *
     *      Alg[x] / poly_modulus
     *
     * @tparam N       Fixed polynomial capacity
     * @tparam PolyMod Polynomial reduction policy (x^N ± 1, none, etc.)
     * @tparam Als     Algebraic structure for coefficient arithmetic (ZmodNField, ZmodNRing)
     * @tparam CoeffType Coefficient type stored in the polynomial
     */
    template<
        std::size_t N,
        typename PolyMod = NoPolyModulus<N>,
        typename CoeffType = OTPQ_ALS_UNDERLYING_TYPE,
        typename Als = arith::ZmodNField<CoeffType> >
        requires core::Field<Als, CoeffType>
    class PolynomialField : public Polynomial<N, CoeffType> {
    public:
        using Base = Polynomial<N, CoeffType>;
        using Base::coeffs;

        Als alg{};
        PolyMod mod{};

        /// Deleted: must provide algebra for initialization.
        PolynomialField() = delete;

        /// Constructs an empty polynomial field with algebra + modulus policy.
        explicit PolynomialField(const Als &a, const PolyMod &m = PolyMod())
            : alg(a), mod(m) {
        }

        /// Constructs and immediately reduces a polynomial.
        explicit PolynomialField(const std::array<CoeffType, N> &arr,
                                 const Als &a = Als(),
                                 const PolyMod &m = PolyMod())
            : Base(arr), alg(a), mod(m) {
            reduce();
        }

        [[nodiscard]] bool is_zero() const {
            return this->degree() == 0 && alg.is_zero(coeffs[0]);
        }


        /// Polynomial addition.
        friend PolynomialField operator+(const PolynomialField &a,
                                         const PolynomialField &b) {
            PolynomialField r(a.alg, a.mod);
            for (std::size_t i = 0; i < N; i++)
                r.coeffs[i] = a.alg.add(a[i], b[i]);
            r.reduce();
            return r;
        }

        /// Polynomial subtraction.
        friend PolynomialField operator-(const PolynomialField &a,
                                         const PolynomialField &b) {
            PolynomialField r(a.alg, a.mod);
            for (std::size_t i = 0; i < N; i++)
                r.coeffs[i] = a.alg.sub(a[i], b[i]);
            r.reduce();
            return r;
        }

        /// Polynomial multiplication (naive convolution with safety check).
        friend PolynomialField operator*(const PolynomialField &a,
                                         const PolynomialField &b) {
            const std::size_t deg_a = a.degree();
            const std::size_t deg_b = b.degree();

            // Step 1: multiply into dynamic buffer (full convolution result)
            std::vector<CoeffType> tmp(deg_a + deg_b + 1, a.alg.zero());

            for (std::size_t i = 0; i <= deg_a; ++i) {
                for (std::size_t j = 0; j <= deg_b; ++j) {
                    tmp[i + j] = a.alg.add(
                        tmp[i + j],
                        a.alg.mul(a[i], b[j])
                    );
                }
            }

            // Step 2: reduce tmp → array<N>
            std::array<CoeffType, N> reduced = a.mod.reduce(tmp, a.alg);

            // Step 3: build final polynomial
            return PolynomialField(reduced, a.alg, a.mod);
        }

    private:
        /// Applies polynomial modulus reduction to this polynomial.
        void reduce() {
            mod.reduce(coeffs, alg);
        }
    };

    /**
     * @brief Cyclotomic polynomial ring Alg[x] / (x^N - 1).
     */
    template<std::size_t N,
        typename CoeffType = OTPQ_ALS_UNDERLYING_TYPE,
        typename Als = arith::ZmodNField<CoeffType> >
    using PolynomialFieldCyclotomic =
    PolynomialField<N, CyclotomicPolyModulus<N>, CoeffType, Als>;
}
