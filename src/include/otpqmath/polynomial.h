#pragma once

#include <array>
#include <cstddef>
#include <vector>


/**
 * @brief Polynomial utilities and polynomial-based ALSs constructions.
 *
 * This namespace defines:
 *   - `Polynomial<N>`: fixed-size polynomial container
 *   - modular polynomial reduction policies
 *       • `NoPolyModulus<N>`               – raw truncation
 *       • `CyclotomicPolyModulus<N>`       – mod (x^N - 1)
 *       • `AntiCyclotomicPolyModulus<N>`   – mod (x^N + 1)
 *   - `PolynomialField<N>`: polynomial ring/field with coefficient
 *     arithmetic delegated to an algebraic structure (e.g., ZmodNField)
 *
 */
namespace otpq::math::als::polynomial {
    /**
     * @brief Fixed-size polynomial with coefficient array of length N.
     *
     * Represents:
     *
     *      c[0] + c[1] x + ... + c[N-1] x^(N-1)
     *
     * The polynomial does not grow dynamically. Reads beyond N return zero;
     * writes beyond N are ignored.
     *
     * @tparam N         Maximum number of coefficients
     * @tparam CoeffType Scalar coefficient type
     */
    template<std::size_t N, typename CoeffType = OTPQ_ALS_UNDERLYING_TYPE>
    class Polynomial {
    public:
        std::array<CoeffType, N> coeffs{}; ///< coefficient storage

        // 1. Default: zero polynomial
        // ------------------------------------------------------------
        constexpr Polynomial() = default;

        // ------------------------------------------------------------
        // 2. Construct from full array of size N
        // ------------------------------------------------------------
        constexpr explicit Polynomial(const std::array<CoeffType, N> &values)
            : coeffs(values) {
        }

        // ------------------------------------------------------------
        // 3. Construct from std::span of size N
        //     (Works with std::array, std::vector, C-arrays, etc.)
        // ------------------------------------------------------------
        constexpr explicit Polynomial(std::span<const CoeffType, N> s) {
            for (std::size_t i = 0; i < N; ++i)
                coeffs[i] = s[i];
        }

        // ------------------------------------------------------------
        // 4. Construct from a container of arbitrary size
        //    Excess values truncated; missing values = 0
        // ------------------------------------------------------------
        template<typename Container>
            requires (!std::is_same_v<std::decay_t<Container>, Polynomial>)
        constexpr explicit Polynomial(const Container &c) {
            std::size_t i = 0;
            for (auto &&v: c) {
                if (i >= N) break;
                coeffs[i++] = static_cast<CoeffType>(v);
            }
            // remaining coeffs default-initialized to zero
        }

        // ------------------------------------------------------------
        // 5. Construct polynomial from initializer_list
        //    Example: Polynomial<5>{1, 2, 3} → [1,2,3,0,0]
        // ------------------------------------------------------------
        constexpr Polynomial(std::initializer_list<CoeffType> list) {
            std::size_t i = 0;
            for (auto &&v: list) {
                if (i >= N) break;
                coeffs[i++] = v;
            }
        }

        // ------------------------------------------------------------
        // 6. Construct constant polynomial:  c(x) = scalar
        // ------------------------------------------------------------
        constexpr explicit Polynomial(CoeffType constant_value) {
            coeffs.fill(CoeffType{});
            coeffs[0] = constant_value;
        }

        /**
         * @brief Read operator with dynamic-style bounds.
         *
         * @return coeff[i] when i < N; otherwise throws exception.
         */
        constexpr CoeffType operator[](std::size_t i) const {
            if (i >= N)
                throw std::out_of_range("[Polynomial] index out of bounds");
            return coeffs[i];
        }

        /**
         * @brief Write operator. Writes beyond N are ignored.
         *
         * @return Reference to coeff[i]; otherwise throws exception.
         */
        constexpr CoeffType &operator[](std::size_t i) {
            if (i >= N)
                throw std::out_of_range("[Polynomial] index out of bounds");
            return coeffs[i];
        }

        /**
         * @brief Returns the static polynomial degree bound (always N).
         *
         * This reflects capacity, not the algebraic degree.
         */
        [[nodiscard]] constexpr std::size_t capacity() const noexcept {
            return N;
        }

        /**
         * @brief Returns the index of the highest non-zero coefficient.
         *
         * If the polynomial is zero, returns 0.
         */
        [[nodiscard]] std::size_t degree() const noexcept {
            for (std::size_t i = N; i-- > 0;) {
                if (coeffs[i] != CoeffType{})
                    return i;
            }
            return 0;
        }
    };

    /**
     * @brief Pretty-prints a polynomial such as:  3 + 5x - 2x^3
     */
    template<std::size_t N, typename CoeffType = OTPQ_ALS_UNDERLYING_TYPE>
    std::ostream &operator<<(std::ostream &os, const Polynomial<N, CoeffType> &p) {
        bool first = true;

        for (std::size_t i = 0; i < N; ++i) {
            CoeffType a = p.coeffs[i];
            if (a == CoeffType{}) continue;

            if (!first) {
                os << (a >= CoeffType{} ? " + " : " - ");
            } else if (a < CoeffType{}) {
                os << "-";
            }

            CoeffType absval = (a < CoeffType{}) ? -a : a;

            if (i == 0 || absval != CoeffType{1})
                os << absval;

            if (i >= 1) {
                os << "x";
                if (i >= 2)
                    os << "^" << i;
            }

            first = false;
        }

        if (first)
            os << "0";

        return os;
    }
}
