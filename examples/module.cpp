#include <iostream>
#include <array>

#include <zmod_arithmetic.h>
#include <polynomial_fields.h>
#include <module.h>

using otpq::math::als::arith::ZmodNField;
using otpq::math::als::polynomial::CyclotomicPolyModulus;
using otpq::math::als::module::PolynomialFieldVector;
using otpq::math::als::module::PolynomialFieldMatrix;

int main() {
    // -----------------------------------------------------------
    // Algebra: R = Z_17[x] / (x^4 - 1)
    // -----------------------------------------------------------
    constexpr int q = 17;
    constexpr std::size_t N = 4;

    ZmodNField alg{q};
    using Mod = CyclotomicPolyModulus<N>;

    // -----------------------------------------------------------
    // Create vector v ∈ R^K from explicit coefficients
    // -----------------------------------------------------------
    constexpr std::size_t K = 2; // rows
    constexpr std::size_t L = 3; // columns

    const std::array<std::array<int, N>, K> vec_coeffs{{
        { {1, 2, 3, 4} },
        { {5, 6, 7, 8} }
    }};

    PolynomialFieldVector v{vec_coeffs, alg, Mod{}};

    std::cout << "Vector v:\n" << v << "\n\n";

    // -----------------------------------------------------------
    // Create matrix A ∈ R^(K × L) from explicit coefficients
    // -----------------------------------------------------------
    const std::array<std::array<std::array<int, N>, L>, K> mat_coeffs{{
        // Row 0
        {{
            { {1,0,0,0} },
            { {0,1,0,0} },
            { {0,0,1,0} }
        }},
        // Row 1
        {{
            { {2,2,2,2} },
            { {3,3,3,3} },
            { {4,4,4,4} }
        }}
    }};

    PolynomialFieldMatrix A{mat_coeffs, alg, Mod{}};

    std::cout << "Matrix A:\n" << A << "\n";

    return 0;
}
