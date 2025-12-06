#include <iostream>
#include <array>
#include <vector>
#include <span>

#include <otpqmath/als/rings/poly/polynomial.h>

using otpq::math::als::poly::polynomial::Polynomial;

int main() {
    std::cout << "=== Polynomial Demonstration ===\n\n";

    // --------------------------------------------------
    // 1. Construct polynomial p(x) = 3 - x + 4x^3  with max degree < 5
    // --------------------------------------------------
    Polynomial<5> p; // zero polynomial: [0,0,0,0,0]

    p[0] = 3; // constant term
    p[1] = -1; // −x
    p[3] = 4; // +4x³

    std::cout << "p(x) = " << p << "\n\n";

    std::cout << "p.capacity()    = " << p.capacity() << "\n";
    std::cout << "p.degree()      = " << p.degree() << "\n\n";


    // --------------------------------------------------
    // 2. Construct q(x) directly from std::array
    //    q(x) = 5 + 2x² + x³
    // --------------------------------------------------
    constexpr std::array init_q{5, 0, 2, 1, 0};
    Polynomial<5> q{init_q};

    std::cout << "q(x) = " << q << "\n\n";


    // --------------------------------------------------
    // 3. Construct r(x) using initializer_list
    //    r(x) = 7 + 3x + 1x²
    // --------------------------------------------------
    Polynomial<5> r{7, 3, 1};

    std::cout << "r(x) = " << r << "\n\n";


    // --------------------------------------------------
    // 4. Construct s(x) from a std::vector
    //    s(x) = 9 + 8x + 7x²
    // --------------------------------------------------
    std::vector<int> vec = {9, 8, 7};
    Polynomial<5> s{vec}; // container constructor

    std::cout << "s(x) = " << s << "\n\n";


    // --------------------------------------------------
    // 5. Construct t(x) as a constant polynomial: t(x) = 42
    // --------------------------------------------------
    Polynomial<5> t{42};

    std::cout << "t(x) = " << t << "\n\n";


    // --------------------------------------------------
    // 6. Demonstrate raw(), unchecked access (e.g., for FFT kernels)
    // --------------------------------------------------
    std::cout << "Using raw() to modify p(x):\n";
    p.raw(2) = -5; // p(x) now has −5x²
    std::cout << "p(x) after raw write: " << p << "\n\n";


    // --------------------------------------------------
    // 7. Demonstrate data(): direct access to underlying storage
    // --------------------------------------------------
    std::cout << "Accessing underlying data() of q:\n";
    auto &qdata = q.data();

    for (std::size_t i = 0; i < q.capacity(); i++)
        std::cout << "  q[" << i << "] = " << qdata[i] << "\n";

    std::cout << "\n";


    // --------------------------------------------------
    // 8. Demonstrate constructing from a span
    // --------------------------------------------------
    std::array<int, 5> span_src{1, 2, 3, 4, 5};
    std::span<const int, 5> span_view{span_src};

    Polynomial<5> u{span_view};

    std::cout << "u(x) constructed from std::span = " << u << "\n\n";


    std::cout << "=== End of Demonstration ===\n";
    return 0;
}
