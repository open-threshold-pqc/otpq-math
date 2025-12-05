#include <iostream>
#include <array>
#include <vector>

#include <otpqmath/numbers/fields/prime_field.h>
#include <otpqmath/poly/fields/polynomial_field.h>

using otpq::math::als::field::PrimeField;
using otpq::math::als::polynomial::PolynomialField;
using otpq::math::als::polynomial::PolynomialFieldCyclotomic;

int main() {
    constexpr std::size_t N = 8;

    using Field = PrimeField<int>;
    using Poly = PolynomialFieldCyclotomic<N>;

    constexpr Field F7{7};

    std::cout << std::format("=== PolynomialFieldCyclotomic (x^{} - 1) Constructor Demo ===\n\n", N);

    // ------------------------------------------------------------
    // 1. Zero polynomial: requires explicit algebra
    // ------------------------------------------------------------
    constexpr Poly a{F7};
    std::cout << "a(x) = " << a << "\n";


    // ------------------------------------------------------------
    // 2. From std::array<N>
    // ------------------------------------------------------------
    constexpr std::array arr_b{1, 2, 3, 4, 5, 6, 0, 0};
    constexpr Poly b{arr_b, F7};
    std::cout << "b(x) = " << b << "\n";


    // ------------------------------------------------------------
    // 3. From span<N>
    // ------------------------------------------------------------
    std::array arr_c{3, 1, 4, 1, 5, 9, 0, 0};
    const std::span<const int, N> sp_c{arr_c};
    const Poly c{sp_c, F7};
    std::cout << "c(x) = " << c << "\n";


    // ------------------------------------------------------------
    // 4. From initializer_list
    //    {1,2,3} → [1,2,3,0,0,0,0,0] (auto padded)
    // ------------------------------------------------------------
    constexpr Poly d{{1, 2, 3}, F7};
    std::cout << "d(x) = " << d << "\n";


    // ------------------------------------------------------------
    // 5. From constant polynomial p(x) = constant_value
    // ------------------------------------------------------------
    constexpr Poly e{42, F7};
    std::cout << "e(x) = " << e << "\n";


    // ------------------------------------------------------------
    // 6. From dynamic container (truncated/padded)
    // ------------------------------------------------------------
    const std::vector<int> vec_f = {9, 8, 7, 6};
    const Poly f{vec_f, F7};
    std::cout << "f(x) = " << f << "\n";

    std::cout << "\n ---------------- \n";
    std::cout << "d(x) + e(x) = " << d - e << "\n";
    std::cout << "b(x) * c(x) = " << b * c << "\n";

    return 0;
}
