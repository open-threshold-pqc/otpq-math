#pragma once

#include <vector>
#include <cstddef>

#include <fields/polynomial_fields.h>

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
        typename PolyMod = polynomial::CyclotomicPolyModulus<N>,
        typename CoeffType = OTPQ_ALS_UNDERLYING_TYPE,
        typename Als = field::PrimeField<CoeffType> >
        requires core::Field<Als, CoeffType>
    class PolynomialFieldVector {
    public:
        using Poly = polynomial::PolynomialField<N, PolyMod, CoeffType, Als>;

        /**
         * @brief Constructs a vector of K zero polynomials.
         *
         * @param alg  Algebraic structure for coefficient arithmetic.
         * @param mod  Polynomial modulus policy.
         */
        explicit PolynomialFieldVector(const Als &alg,
                                       const PolyMod &mod = PolyMod()) {
            data_.reserve(K);
            for (std::size_t i = 0; i < K; ++i)
                data_.emplace_back(alg, mod);
        }

        /**
         * @brief Constructs from K arrays of N coefficients.
         *
         * @param coeffs  Array of K arrays of N coefficients.
         * @param alg     Algebraic structure.
         * @param mod     Polynomial modulus policy.
         */
        explicit PolynomialFieldVector(
            const std::array<std::array<CoeffType, N>, K> &coeffs,
            const Als &alg,
            const PolyMod &mod = PolyMod()) {
            data_.reserve(K);
            for (std::size_t i = 0; i < K; ++i)
                data_.emplace_back(coeffs[i], alg, mod);
        }

        /**
         * @brief Constructs from a span of K arrays of N coefficients.
         *
         * @param spans   Span of K coefficient arrays.
         * @param alg     Algebraic structure.
         * @param mod     Polynomial modulus policy.
         */
        explicit PolynomialFieldVector(
            std::span<const std::array<CoeffType, N>, K> spans,
            const Als &alg,
            const PolyMod &mod = PolyMod()) {
            data_.reserve(K);
            for (std::size_t i = 0; i < K; ++i)
                data_.emplace_back(spans[i], alg, mod);
        }

        /**
         * @brief Constructs from an arbitrary container of K coefficient-containers.
         *
         * Each element of @c container[i] may have any size; elements are truncated
         * or padded to length N.
         *
         * @tparam Container  Any container with K sub-containers yielding coefficients.
         * @param container   Nested coefficient container.
         * @param alg         Algebraic structure.
         * @param mod         Polynomial modulus policy.
         */
        template<typename Container>
            requires (!std::is_same_v<std::decay_t<Container>, PolynomialFieldVector>)
        explicit PolynomialFieldVector(const Container &container,
                                       const Als &alg,
                                       const PolyMod &mod = PolyMod()) {
            data_.reserve(K);
            std::size_t i = 0;

            for (auto &&c: container) {
                if (i >= K) break;
                data_.emplace_back(c, alg, mod); // uses container-based Poly constructor
                ++i;
            }

            for (; i < K; ++i)
                data_.emplace_back(alg, mod);
        }

        /**
         * @brief Constructs each component as a constant polynomial.
         *
         * @param constant  Constant coefficient for each component polynomial.
         * @param alg       Algebraic structure.
         * @param mod       Polynomial modulus policy.
         */
        explicit PolynomialFieldVector(CoeffType constant,
                                       const Als &alg,
                                       const PolyMod &mod = PolyMod()) {
            data_.reserve(K);
            for (std::size_t i = 0; i < K; ++i)
                data_.emplace_back(constant, alg, mod);
        }

        /**
         * @brief Constructs from initializer_list of polynomials.
         *
         * Fewer than K → remaining are zero polynomials.
         * More than K → excess elements ignored.
         *
         * @param list  List of PolynomialField objects.
         * @param alg   Algebraic structure (must match input polynomials).
         * @param mod   Polynomial modulus policy.
         */
        explicit PolynomialFieldVector(std::initializer_list<Poly> list,
                                       const Als &alg,
                                       const PolyMod &mod = PolyMod()) {
            data_.reserve(K);
            std::size_t i = 0;

            for (auto &&p: list) {
                if (i >= K) break;
                data_.emplace_back(p);
                ++i;
            }

            for (; i < K; ++i)
                data_.emplace_back(alg, mod);
        }

        /**
         * @brief Returns the module dimension K.
         */
        [[nodiscard]] constexpr std::size_t size() const noexcept { return K; }

        /**
         * @brief Read/write access to the i-th polynomial component.
         */
        Poly &operator[](std::size_t i) { return data_[i]; }

        /**
         * @brief Read-only access to the i-th polynomial component.
         */
        const Poly &operator[](std::size_t i) const { return data_[i]; }

        /**
         * @brief Elementwise vector addition.
         */
        friend PolynomialFieldVector operator+(
            const PolynomialFieldVector &a,
            const PolynomialFieldVector &b) {
            PolynomialFieldVector r(a.data_[0].algebra(), a.data_[0].modulus());

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
            PolynomialFieldVector r(a.data_[0].algebra(), a.data_[0].modulus());

            for (std::size_t i = 0; i < K; ++i)
                r[i] = a[i] - b[i];

            return r;
        }

    private:
        /**
         * @brief Underlying storage for the K component polynomials.
         *
         * std::vector is used because PolynomialField has no default constructor,
         * and each element must be explicitly constructed with (Alg, Mod).
         */
        std::vector<Poly> data_;
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
     * @brief Left scalar multiplication by a polynomial.
     *
     * Computes the product
     *      r[i] = a * v[i]
     * for each component of the vector, producing an element of R^K where
     * R = Als[x] / PolyMod.  Both operands must share the same algebra and
     * modulus policy.
     *
     * @tparam K         Vector length
     * @tparam N         Polynomial capacity
     * @tparam PolyMod   Polynomial modulus policy
     * @tparam CoeffType Coefficient type
     * @tparam Als       Coefficient algebra
     *
     * @param a  Scalar polynomial in R
     * @param v  Vector over R
     *
     * @return Vector whose entries are a · v[i]
     */
    template<std::size_t K,
        std::size_t N,
        typename PolyMod,
        typename CoeffType,
        typename Als>
    PolynomialFieldVector<K, N, PolyMod, CoeffType, Als>
    operator*(const polynomial::PolynomialField<N, PolyMod, CoeffType, Als> &a,
              const PolynomialFieldVector<K, N, PolyMod, CoeffType, Als> &v) {
        PolynomialFieldVector<K, N, PolyMod, CoeffType, Als> r(a.algebra(), a.modulus());
        for (std::size_t i = 0; i < K; ++i)
            r[i] = a * v[i];
        return r;
    }

    /**
     * @brief Right scalar multiplication by a polynomial.
     *
     * Computes the product
     *      r[i] = v[i] * a
     * for each component of the vector.  Equivalent to left multiplication
     * when multiplication in R is commutative.
     *
     * @tparam K         Vector length
     * @tparam N         Polynomial capacity
     * @tparam PolyMod   Polynomial modulus policy
     * @tparam CoeffType Coefficient type
     * @tparam Als       Coefficient algebra
     *
     * @param v  Vector over R
     * @param a  Scalar polynomial in R
     *
     * @return Vector whose entries are v[i] · a
     */
    template<std::size_t K,
        std::size_t N,
        typename PolyMod,
        typename CoeffType,
        typename Als>
    PolynomialFieldVector<K, N, PolyMod, CoeffType, Als>
    operator*(const PolynomialFieldVector<K, N, PolyMod, CoeffType, Als> &v,
              const polynomial::PolynomialField<N, PolyMod, CoeffType, Als> &a) {
        PolynomialFieldVector<K, N, PolyMod, CoeffType, Als> r(a.algebra(), a.modulus());
        for (std::size_t i = 0; i < K; ++i)
            r[i] = v[i] * a;
        return r;
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
        typename PolyMod = polynomial::CyclotomicPolyModulus<N>,
        typename CoeffType = OTPQ_ALS_UNDERLYING_TYPE,
        typename Als = field::PrimeField<CoeffType> >
        requires core::Field<Als, CoeffType>
    class PolynomialFieldMatrix {
    public:
        using Poly = polynomial::PolynomialField<N, PolyMod, CoeffType, Als>;
        using VecCols = PolynomialFieldVector<L, N, PolyMod, CoeffType, Als>;
        using VecRows = PolynomialFieldVector<K, N, PolyMod, CoeffType, Als>;

        /**
     * @brief Constructs a K×L zero matrix.
     */
        explicit PolynomialFieldMatrix(const Als &alg,
                                       const PolyMod &mod = PolyMod()) {
            data_.resize(K);
            for (std::size_t i = 0; i < K; ++i) {
                data_[i].reserve(L);
                for (std::size_t j = 0; j < L; ++j)
                    data_[i].emplace_back(alg, mod);
            }
        }

        /**
         * @brief Constructs from K×L arrays of N coefficients.
         */
        explicit PolynomialFieldMatrix(
            const std::array<std::array<std::array<CoeffType, N>, L>, K> &coeffs,
            const Als &alg,
            const PolyMod &mod = PolyMod()) {
            data_.resize(K);
            for (std::size_t i = 0; i < K; ++i) {
                data_[i].reserve(L);
                for (std::size_t j = 0; j < L; ++j)
                    data_[i].emplace_back(coeffs[i][j], alg, mod);
            }
        }

        /**
         * @brief Constructs from a span of K rows of L coefficient arrays.
         */
        explicit PolynomialFieldMatrix(
            std::span<const std::array<std::array<CoeffType, N>, L>, K> rows,
            const Als &alg,
            const PolyMod &mod = PolyMod()) {
            data_.resize(K);
            for (std::size_t i = 0; i < K; ++i) {
                data_[i].reserve(L);
                for (std::size_t j = 0; j < L; ++j)
                    data_[i].emplace_back(rows[i][j], alg, mod);
            }
        }

        /**
         * @brief Constructs from an arbitrary nested container.
         */
        template<class Container>
            requires (!std::is_same_v<std::decay_t<Container>, PolynomialFieldMatrix>)
        explicit PolynomialFieldMatrix(const Container &rows,
                                       const Als &alg,
                                       const PolyMod &mod = PolyMod()) {
            data_.resize(K);
            std::size_t i = 0;

            for (auto &&row: rows) {
                if (i >= K) break;

                data_[i].reserve(L);
                std::size_t j = 0;

                for (auto &&cell: row) {
                    if (j >= L) break;
                    data_[i].emplace_back(cell, alg, mod);
                    ++j;
                }

                for (; j < L; ++j)
                    data_[i].emplace_back(alg, mod);

                ++i;
            }

            for (; i < K; ++i) {
                data_[i].reserve(L);
                for (std::size_t j = 0; j < L; ++j)
                    data_[i].emplace_back(alg, mod);
            }
        }

        /**
         * @brief Constructs all entries as constant polynomials.
         */
        explicit PolynomialFieldMatrix(CoeffType constant,
                                       const Als &alg,
                                       const PolyMod &mod = PolyMod()) {
            data_.resize(K);
            for (std::size_t i = 0; i < K; ++i) {
                data_[i].reserve(L);
                for (std::size_t j = 0; j < L; ++j)
                    data_[i].emplace_back(constant, alg, mod);
            }
        }


        /**
         * @brief Read/write access to the [i,j] entry of the matrix.
         */
        Poly &operator[](std::size_t i, std::size_t j) { return data_[i][j]; }

        /**
         * @brief Read-only access to the [i,j] entry of the matrix.
         */
        const Poly &operator[](std::size_t i, std::size_t j) const { return data_[i][j]; }

        /**
         * @brief Matrix addition.
         *
         * Returns A + B computed elementwise.
         */
        friend PolynomialFieldMatrix operator+(
            const PolynomialFieldMatrix &A,
            const PolynomialFieldMatrix &B) {
            PolynomialFieldMatrix R(A.data_[0][0].algebra(), A.data_[0][0].modulus());

            for (std::size_t i = 0; i < K; ++i)
                for (std::size_t j = 0; j < L; ++j)
                    R[i, j] = A[i, j] + B[i, j];

            return R;
        }


        /**
         * @brief Matrix subtraction.
         *
         * Returns A - B computed elementwise.
         */
        friend PolynomialFieldMatrix operator-(
            const PolynomialFieldMatrix &A,
            const PolynomialFieldMatrix &B) {
            PolynomialFieldMatrix R(A.data_[0][0].algebra(), A.data_[0][0].modulus());

            for (std::size_t i = 0; i < K; ++i)
                for (std::size_t j = 0; j < L; ++j)
                    R[i, j] = A[i, j] - B[i, j];

            return R;
        }


        /**
         * @brief Matrix–vector multiplication.
         *
         * A is K×L, x is length L, result is length K.
         */
        friend VecRows operator*(const PolynomialFieldMatrix &A,
                                 const VecCols &x) {
            VecRows y{A[0,0].algebra(), A[0,0].modulus()}; // zero vector of length K

            for (std::size_t i = 0; i < K; ++i) {
                for (std::size_t j = 0; j < L; ++j)
                    y[i] = y[i] + (A[i, j] * x[j]);
            }

            return y;
        }

    private:
        /**
         * @brief Row-major storage of size K×L.
         *
         * Each row is a std::vector<Poly> allowing explicit construction
         * without a default constructor for Poly.
         */
        std::vector<std::vector<Poly> > data_;
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
