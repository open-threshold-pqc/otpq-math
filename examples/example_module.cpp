#include <iostream>
#include <array>

#include <otpqmath/numbers/zmod_arithmetic.h>
#include <otpqmath/poly/module.h>

using otpq::math::als::arith::ZmodNField;
using otpq::math::als::module::PolynomialFieldVector;
using otpq::math::als::module::PolynomialFieldMatrix;

int main() {
    constexpr std::size_t N = 4;
    constexpr std::size_t K = 2;
    constexpr std::size_t L = 3;

    constexpr ZmodNField<> F17{17};

    // -----------------------------------------------------------
    // 1. Create vector v ∈ R^K from explicit coefficients
    // -----------------------------------------------------------
    constexpr std::array<std::array<int, N>, K> vec_coeffs{
        {
            {{1, 2, 3, 4}},
            {{5, 6, 7, 8}}
        }
    };
    const PolynomialFieldVector v{vec_coeffs, F17};
    std::cout << "Vector v:\n" << v << "\n\n";

    // -----------------------------------------------------------
    // 2. Construct another vector using a constant polynomial
    // -----------------------------------------------------------
    const PolynomialFieldVector<L, N> w{3, F17};

    std::cout << "Vector w:\n" << w << "\n\n";


    // -----------------------------------------------------------
    // 3. Construct vector from arbitrary container type
    // -----------------------------------------------------------
    const std::vector<std::vector<int> > coeff_vec = {
        {9, 8, 7, 6},
        {1, 1, 2, 3}
    };

    PolynomialFieldVector<K, N> u{coeff_vec, F17};


    std::cout << "Vector u:\n" << w - w << "\n\n";


    using Mat = PolynomialFieldMatrix<K, L, N>;
    // -----------------------------------------------------------
    // 1. Construct matrix A from explicit coefficient arrays
    // -----------------------------------------------------------
    constexpr std::array<std::array<std::array<int, N>, L>, K> coeff_A{
        {
            {
                {
                    {{1, 2, 3, 4}},
                    {{5, 0, 1, 0}},
                    {{7, 7, 7, 7}}
                }
            },
            {
                {
                    {{2, 4, 6, 8}},
                    {{3, 1, 4, 1}},
                    {{0, 0, 0, 0}}
                }
            }
        }
    };

    Mat A{coeff_A, F17};

    std::cout << "Matrix A:\n" << A << "\n\n";

    // -----------------------------------------------------------
    // 2. Construct B using constant polynomials
    // -----------------------------------------------------------
    const Mat B{3, F17};

    std::cout << "Matrix B:\n" << B << "\n\n";

    // -----------------------------------------------------------
    // 3. Construct matrix C from arbitrary nested container
    // -----------------------------------------------------------
    const std::vector<std::vector<std::vector<int> > > coeff_C = {
        {
            {9, 8, 7, 6}, // row 0, col 0
            {1, 1, 2, 3}, // row 0, col 1
            {5, 5, 5, 5} // row 0, col 2
        },
        {
            {4, 3, 2, 1}, // row 1, col 0
            {6, 7, 8, 9}, // row 1, col 1
            {0} // row 1, col 2 (auto padded)
        }
    };

    const Mat C{coeff_C, F17};

    std::cout << "Matrix C:\n" << C << "\n\n";

    std::cout << "\n -------- \n";
    std::cout << "w : " << w << "\n";
    std::cout << "A : " << A << "\n";
    std::cout << "Aw : " << A * w << "\n";

    return 0;
}
