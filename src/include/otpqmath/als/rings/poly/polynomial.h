#pragma once

#include <array>
#include <cstddef>
#include <span>


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
namespace otpq::math::als::rings::polynomial {
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
    template<std::size_t N, std::integral CoeffType = OTPQ_ALS_UNDERLYING_TYPE>
    class Polynomial {
    public:
        /**
         * @brief Constructs the zero polynomial.
         */
        constexpr Polynomial() = default;

        /**
         * @brief Constructs from an array of exactly N coefficients.
         *
         * @param values  The full coefficient array.
         */
        constexpr explicit Polynomial(const std::array<CoeffType, N> &values)
            : coeffs_(values) {
        }

        /**
         * @brief Constructs from a fixed-size span of exactly N coefficients.
         *
         * @param s  Span containing N coefficients.
         */
        constexpr explicit Polynomial(std::span<const CoeffType, N> s) {
            for (std::size_t i = 0; i < N; ++i)
                coeffs_[i] = s[i];
        }

        /**
         * @brief Constructs from an arbitrary container.
         *
         * Elements are copied into the polynomial:
         *   - up to N elements are stored,
         *   - excess elements are ignored,
         *   - missing elements are initialized to zero.
         *
         * This constructor supports all input ranges (vectors, arrays, lists, etc.)
         * except another Polynomial.
         *
         * @param c  Input container or range.
         */
        template<typename Container>
            requires (!std::is_same_v<std::decay_t<Container>, Polynomial>)
        constexpr explicit Polynomial(const Container &c) {
            std::size_t i = 0;
            for (auto &&v: c) {
                if (i >= N) break;
                coeffs_[i++] = static_cast<CoeffType>(v);
            }
        }

        /**
         * @brief Constructs from an initializer list.
         *
         * Example: `Polynomial<5>{1,2,3}` becomes `[1,2,3,0,0]`.
         *
         * @param list  Coefficient list.
         */
        constexpr Polynomial(std::initializer_list<CoeffType> list) {
            std::size_t i = 0;
            for (auto &&v: list) {
                if (i >= N) break;
                coeffs_[i++] = v;
            }
        }

        /**
         * @brief Constructs a constant polynomial p(x) = constant_value.
         *
         * @param constant_value Constant coefficient.
         */
        constexpr explicit Polynomial(CoeffType constant_value) {
            coeffs_.fill(CoeffType{});
            coeffs_[0] = constant_value;
        }

        /**
     * @brief Pretty-prints a polynomial such as: 3 + 5x - 2x^3
     */
        friend std::ostream &operator<<(std::ostream &os,
                                        const Polynomial &p) {
            bool first = true;

            for (std::size_t i = 0; i < N; ++i) {
                CoeffType a = p.coeffs_[i];
                if (a == CoeffType{}) continue;

                if (!first)
                    os << (a >= CoeffType{} ? " + " : " - ");
                else if (a < CoeffType{})
                    os << "-";

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

        /**
         * @brief Read-only coefficient access with bounds checking.
         *
         * @param i  Coefficient index.
         * @return   coeff[i]
         *
         * @throws std::out_of_range if i ≥ N.
         */
        constexpr const CoeffType &operator[](std::size_t i) const {
            if (i >= N)
                throw std::out_of_range("[Polynomial] index out of bounds");
            return coeffs_[i];
        }

        /**
         * @brief Writable coefficient access with bounds checking.
         *
         * @param i  Coefficient index.
         * @return   Reference to coeff[i]
         *
         * @throws std::out_of_range if i ≥ N.
         */
        [[nodiscard]]
        constexpr CoeffType &operator[](std::size_t i) {
            if (i >= N)
                throw std::out_of_range("[Polynomial] index out of bounds");
            return coeffs_[i];
        }

        /**
         * @brief Returns a direct, unchecked reference to the i-th coefficient.
         *
         * This function performs **no bounds checking** and is intended for
         * performance-critical routines such as NTT/FFT kernels or SIMD loops
         * where index validity is guaranteed externally.
         *
         * @param i  Coefficient index (must satisfy 0 ≤ i < N).
         * @return   Reference to coeffs_[i] (unchecked).
         */
        constexpr CoeffType &raw(std::size_t i) noexcept { return coeffs_[i]; }

        /**
         * @brief Returns a direct, unchecked reference to the i-th coefficient
         *        (const overload).
         *
         * See the non-const overload for details.
         */
        constexpr const CoeffType &raw(std::size_t i) const noexcept { return coeffs_[i]; }


        /**
         * @brief Provides direct access to the underlying coefficient array.
         *
         * Useful for:
         *   • passing coefficients to NTT/FFT routines,
         *   • SIMD/vectorized operations,
         *   • interoperability with low-level C APIs.
         *
         * This returns the internal std::array<N> by reference; no copy occurs.
         *
         * @return Reference to the underlying coefficient storage.
         */
        constexpr auto &data() noexcept { return coeffs_; }

        /**
         * @brief Const-qualified access to the underlying coefficient array.
         *
         * See the non-const overload for details.
         */
        constexpr const auto &data() const noexcept { return coeffs_; }

        /**
         * @brief Returns the static coefficient capacity N.
         *
         * This reflects storage capacity, not algebraic degree.
         *
         * @return N.
         */
        [[nodiscard]] constexpr std::size_t capacity() const noexcept { return N; }

        /**
         * @brief Returns the polynomial's algebraic degree.
         *
         * The degree is the highest index i such that cᵢ ≠ 0.
         * If all coefficients are zero, the degree is defined as 0.
         *
         * @return Degree index, or 0 if polynomial is identically zero.
         */
        [[nodiscard]] std::size_t degree() const noexcept {
            for (std::size_t i = N; i-- > 0;)
                if (coeffs_[i] != CoeffType{})
                    return i;
            return 0;
        }

    protected:
        /**
         * @brief Protected coefficient storage.
         *
         * Made protected so derived classes (e.g. PolynomialField) can manipulate
         * coefficients exactly as before, without API changes, while keeping
         * encapsulation for outside consumers.
         */
        std::array<CoeffType, N> coeffs_{};
    };
}
