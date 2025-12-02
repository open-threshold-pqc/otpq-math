#include <iostream>
#include <array>
#include "polynomial.h"

using otpq::math::als::polynomial::Polynomial;

int main() {
    std::cout << "=== Polynomial Demonstration ===\n\n";

    // --------------------------------------------------
    // Construct polynomial p(x) = 3 - x + 4x^3  with max degree < 5
    // --------------------------------------------------
    Polynomial<5> p; // all coefficients begin as zero

    p[0] = 3; // constant term
    p[1] = -1; // −x
    p[3] = 4; // +4x³

    std::cout << "p(x) = " << p << "\n\n";

    std::cout << "p capacity    = " << p.capacity() << "\n";
    std::cout << "actual degree = " << p.degree() << "\n\n";

    // --------------------------------------------------
    // Construct polynomial q(x) from an array
    // q(x) = 5 + 2x² + x³
    // --------------------------------------------------
    constexpr std::array init_q{5, 0, 2, 1, 0};
    Polynomial<5> q{init_q};

    std::cout << "q(x) = " << q << "\n\n";

    std::cout << "q capacity    = " << q.capacity() << "\n";
    std::cout << "actual degree = " << q.degree() << "\n\n";

    // --------------------------------------------------
    // Print coefficients explicitly
    // --------------------------------------------------
    std::cout << "q(x) coefficients:\n";
    for (std::size_t i = 0; i <= q.degree(); ++i)
        std::cout << "  c[" << i << "] = " << q[i] << "\n";


    std::cout << "\nDone.\n";


    return 0;
}
