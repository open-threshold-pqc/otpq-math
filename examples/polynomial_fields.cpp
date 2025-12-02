#include <iostream>
#include <array>

#include <zmod_arithmetic.h>
#include <polynomial_fields.h>

int main() {
    using otpq::math::als::arith::ZmodNField;
    using otpq::math::als::polynomial::PolynomialField;
    using otpq::math::als::polynomial::PolynomialFieldCyclotomic;

    std::cout << "=== PolynomialField (No Modulus) Example ===\n\n";

    constexpr std::size_t N = 10;

    // Coefficient field F₁₁
    const ZmodNField F11{11};
    std::cout << "Coefficient field modulus: " << F11.mod << "\n\n";

    using Poly = PolynomialField<N>;

    // -------------------------
    // p(x) = 3 + 5x + 7x³
    // -------------------------
    Poly p{
        std::array{3, 5, 0, 7, 0, 0, 0, 0, 0, 0},
        F11
    };

    // -------------------------
    // q(x) = 4 + x² + 9x⁵
    // -------------------------
    Poly q{
        std::array{4, 0, 1, 0, 0, 9, 0, 0, 0, 0},
        F11
    };

    std::cout << "p(x) = " << p << "\n";
    std::cout << "q(x) = " << q << "\n\n";

    std::cout << "p + q = " << (p + q) << "\n";
    std::cout << "p - q = " << (p - q) << "\n";

    // Multiplication with safety check
    try {
        auto r = p * q;
        std::cout << "p * q = " << r << "\n\n";
    } catch (const std::exception& e) {
        std::cout << "Multiplication error: " << e.what() << "\n\n";
    }

    // ===========================================================
    // Cyclotomic Example (mod x^N - 1)
    // ===========================================================

    std::cout << "=== PolynomialFieldCyclotomic Example (x^N - 1) ===\n\n";

    using PolyCyc = PolynomialFieldCyclotomic<N>;

    const PolyCyc f{
        std::array{4, 0, 1, 0, 0, 9, 0, 0, 0, 0},
        11
    };

    const PolyCyc g{
        std::array{4, 0, 1, 0, 0, 9, 0, 0, 0, 0},
        11
    };

    std::cout << "f(x) = " << f << "\n";
    std::cout << "g(x) = " << g << "\n";
    std::cout << "f * g (x) = " << (f * g) << "\n";


    return 0;
}
