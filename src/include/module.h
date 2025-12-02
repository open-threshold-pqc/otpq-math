#pragma once

#include <vector>
#include <cstddef>

#include <polynomial_fields.h>

namespace otpq::math::als::module {

    /**
     * @brief A fixed-size module vector over a polynomial ring/field.
     *
     * This class represents an element of the free module
     *    R^K   where   R = Als[x] / PolyMod
     *
     * Each component is a PolynomialField<N>, and the vector supports:
     *   • elementwise vector addition and subtraction
     *   • direct access via operator[]
     *   • construction from (alg, mod) or explicit coefficient arrays
     *
     * @tparam K         Number of polynomial components (vector dimension)
     * @tparam N         Degree bound of each polynomial (size of coefficient array)
     * @tparam PolyMod   Modulus policy for polynomial reduction (e.g., x^N − 1)
     * @tparam CoeffType Underlying scalar coefficient type (e.g., int)
     * @tparam Als       Algebraic structure implementing addition/multiplication
     *                   over CoeffType (e.g., ZmodNField)
     */
    template<
        std::size_t K,
        std::size_t N,
        typename PolyMod = polynomial::NoPolyModulus<N>,
        typename CoeffType = OTPQ_ALS_UNDERLYING_TYPE,
        typename Als = arith::ZmodNField<CoeffType> >
        requires core::Field<Als, CoeffType>
    class PolynomialFieldVector {
    public:
        using Poly = polynomial::PolynomialField<N, PolyMod, CoeffType, Als>;

        /**
         * @brief Storage for the K component polynomials.
         *
         * std::vector is used instead of std::array because PolynomialField
         * does not have a default constructor, and vector::emplace_back()
         * allows explicit construction of each element.
         */
        std::vector<Poly> data;

        /**
         * @brief Constructs a vector of K zero polynomials.
         *
         * Each polynomial is initialized using the provided algebraic structure
         * and modulus policy. No coefficients are set explicitly (all zero).
         *
         * @param alg  Algebraic structure for coefficient arithmetic.
         * @param mod  Polynomial modulus policy.
         */
        PolynomialFieldVector(const Als &alg = Als(),
                              const PolyMod &mod = PolyMod()) {
            data.reserve(K);
            for (std::size_t i = 0; i < K; ++i)
                data.emplace_back(alg, mod);
        }

        /**
         * @brief Constructs a vector from explicit coefficient arrays.
         *
         * Each vector component is initialized from an array of N coefficients,
         * which is reduced according to PolyMod.
         *
         * @param coeffs  Array of K arrays of N coefficients.
         * @param alg     Algebraic structure for coefficient arithmetic.
         * @param mod     Polynomial modulus policy.
         */
        explicit PolynomialFieldVector(
            const std::array<std::array<CoeffType, N>, K> &coeffs,
            const Als &alg = Als(),
            const PolyMod &mod = PolyMod()) {
            data.reserve(K);
            for (std::size_t i = 0; i < K; ++i)
                data.emplace_back(coeffs[i], alg, mod);
        }

        /**
         * @brief Returns the module dimension K.
         */
        [[nodiscard]] constexpr std::size_t size() const noexcept { return K; }

        /**
         * @brief Read/write access to the i-th polynomial component.
         */
        Poly &operator[](std::size_t i) { return data[i]; }

        /**
         * @brief Read-only access to the i-th polynomial component.
         */
        const Poly &operator[](std::size_t i) const { return data[i]; }

        /**
         * @brief Elementwise vector addition.
         */
        friend PolynomialFieldVector operator+(
            const PolynomialFieldVector &a,
            const PolynomialFieldVector &b) {

            PolynomialFieldVector r(a.data[0].alg, a.data[0].mod);

            for (std::size_t i = 0; i < K; ++i)
                r[i] = a[i] + b[i];

            return r;
        }

        /**
         * @brief Elementwise vector subtraction.
         */
        friend PolynomialFieldVector operator-(
            const PolynomialFieldVector &a,
            const PolynomialFieldVector &b) {

            PolynomialFieldVector r(a.data[0].alg, a.data[0].mod);

            for (std::size_t i = 0; i < K; ++i)
                r[i] = a[i] - b[i];

            return r;
        }
    };


    /**
     * @brief Pretty-printer for PolynomialFieldVector.
     *
     * Prints the vector in the form:
     *
     *    [
     *      p0 ,
     *      p1 ,
     *      ... ,
     *      p(K-1)
     *    ]
     */
    template<
        std::size_t K,
        std::size_t N,
        typename PolyMod,
        typename CoeffType,
        typename Als>
    std::ostream &operator<<(
        std::ostream &os,
        const PolynomialFieldVector<K, N, PolyMod, CoeffType, Als> &v) {

        os << "[ \n";
        for (std::size_t i = 0; i < K; ++i) {
            os << "  " << v[i];
            if (i + 1 < K) os << " ,\n";
        }
        os << "\n]";
        return os;
    }


    /**
     * @brief A K×L matrix over the polynomial ring R = Als[x] / PolyMod.
     *
     * This represents a module homomorphism:
     *       R^L  →  R^K
     *
     * Internally stored as a vector of K rows, each containing L polynomials.
     * Supports:
     *   • matrix-vector multiplication
     *   • matrix addition
     *   • element access via operator(i,j)
     *
     * @tparam K         Number of rows
     * @tparam L         Number of columns
     * @tparam N         Degree bound for each polynomial coefficient array
     * @tparam PolyMod   Polynomial modulus policy
     * @tparam CoeffType Scalar coefficient type
     * @tparam Als       Coefficient algebra (e.g., ZmodNField)
     */
    template<
        std::size_t K,
        std::size_t L,
        std::size_t N,
        typename PolyMod = polynomial::NoPolyModulus<N>,
        typename CoeffType = OTPQ_ALS_UNDERLYING_TYPE,
        typename Als = arith::ZmodNField<CoeffType> >
        requires core::Field<Als, CoeffType>
    class PolynomialFieldMatrix {
    public:
        using Poly = polynomial::PolynomialField<N, PolyMod, CoeffType, Als>;
        using VecCols = PolynomialFieldVector<L, N, PolyMod, CoeffType, Als>;
        using VecRows = PolynomialFieldVector<K, N, PolyMod, CoeffType, Als>;

        /**
         * @brief Row-major storage of size K×L.
         *
         * Each row is a std::vector<Poly> allowing explicit construction
         * without a default constructor for Poly.
         */
        std::vector<std::vector<Poly> > data;

        /**
         * @brief Constructs a K×L zero matrix.
         *
         * Each entry is initialized as a zero PolynomialField created with
         * the provided algebraic structure and modulus policy.
         */
        PolynomialFieldMatrix(const Als &alg = Als(),
                              const PolyMod &mod = PolyMod()) {

            data.resize(K);
            for (std::size_t i = 0; i < K; ++i) {
                data[i].reserve(L);
                for (std::size_t j = 0; j < L; ++j)
                    data[i].emplace_back(alg, mod);
            }
        }

        /**
         * @brief Constructs a matrix from coefficient arrays.
         *
         * @param coeffs  A K×L array of coefficient arrays of length N.
         * @param alg     Algebraic structure for coefficient arithmetic.
         * @param mod     Polynomial modulus policy.
         */
        explicit PolynomialFieldMatrix(
            const std::array<std::array<std::array<CoeffType, N>, L>, K> &coeffs,
            const Als &alg = Als(),
            const PolyMod &mod = PolyMod()) {

            data.resize(K);
            for (std::size_t i = 0; i < K; ++i) {
                data[i].reserve(L);
                for (std::size_t j = 0; j < L; ++j)
                    data[i].emplace_back(coeffs[i][j], alg, mod);
            }
        }

        /**
         * @brief Read/write access to the [i,j] entry of the matrix.
         */
        Poly &operator[](std::size_t i, std::size_t j) { return data[i][j]; }

        /**
         * @brief Read-only access to the [i,j] entry of the matrix.
         */
        const Poly &operator[](std::size_t i, std::size_t j) const { return data[i][j]; }

        /**
         * @brief Matrix addition.
         *
         * Returns A + B computed elementwise.
         */
        friend PolynomialFieldMatrix operator+(
            const PolynomialFieldMatrix &A,
            const PolynomialFieldMatrix &B) {

            PolynomialFieldMatrix R(A.data[0][0].alg, A.data[0][0].mod);

            for (std::size_t i = 0; i < K; ++i)
                for (std::size_t j = 0; j < L; ++j)
                    R(i, j) = A(i, j) + B(i, j);

            return R;
        }
    };


    /**
     * @brief Pretty-printer for PolynomialFieldMatrix.
     *
     * Printed format:
     *
     * [
     *   [ a00 , a01 , ... , a0(L-1) ],
     *   [ a10 , a11 , ... , a1(L-1) ],
     *   ...
     *   [ a(K-1)0 , ... , a(K-1)(L-1) ]
     * ]
     */
    template<
        std::size_t K,
        std::size_t L,
        std::size_t N,
        typename PolyMod,
        typename CoeffType,
        typename Als>
    std::ostream &operator<<(
        std::ostream &os,
        const PolynomialFieldMatrix<K, L, N, PolyMod, CoeffType, Als> &M) {

        os << "[\n";
        for (std::size_t i = 0; i < K; ++i) {
            os << "  [ ";
            for (std::size_t j = 0; j < L; ++j) {
                os << M[i, j];
                if (j + 1 < L) os << " , ";
            }
            os << " ]";
            if (i + 1 < K) os << ",";
            os << "\n";
        }
        os << "]";
        return os;
    }

}