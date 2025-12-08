#include <iostream>
#include <array>
#include <vector>
#include <format>

#include <otpqmath/als/field/numbers/prime_field.h>
#include <otpqmath/als/ring/poly/polynomial_ring.h>

using otpq::math::als::field::number::PrimeField;
using otpq::math::als::ring::polynomial::PolynomialRing;
using otpq::math::als::ring::polynomial::PolynomialCyclotomicRing;

int main() {
    constexpr std::size_t N = 8;

    using Field = PrimeField<int>;

    using Poly = PolynomialCyclotomicRing<N>;      // = PolynomialRing<N, PolyCyclotomicModulus<N>>

    constexpr Field F7{7};

    std::cout << std::format("=== PolynomialCyclotomicRing (x^{} - 1) Constructor Demo ===\n\n", N);

    // ------------------------------------------------------------
    // 1. Zero polynomial: requires explicit algebra
    // ------------------------------------------------------------
    constexpr Poly a{F7};
    std::cout << "a(x) = " << a << "\n";

    // ------------------------------------------------------------
    // 2. From std::array<N>
    // ------------------------------------------------------------
    constexpr std::array<int, N> arr_b{1, 2, 3, 4, 5, 6, 0, 0};
    constexpr Poly b{arr_b, F7};
    std::cout << "b(x) = " << b << "\n";

    // ------------------------------------------------------------
    // 3. From span<N>
    // ------------------------------------------------------------
    std::array<int, N> arr_c{3, 1, 4, 1, 5, 9, 0, 0};
    const std::span<const int, N> sp_c{arr_c};
    const Poly c{sp_c, F7};
    std::cout << "c(x) = " << c << "\n";

    // ------------------------------------------------------------
    // 4. From initializer_list {1,2,3}
    //    → padded to [1,2,3,0,0,0,0,0]
    // ------------------------------------------------------------
    constexpr Poly d{{1, 2, 3}, F7};
    std::cout << "d(x) = " << d << "\n";

    // ------------------------------------------------------------
    // 5. Constant polynomial p(x) = constant_value
    // ------------------------------------------------------------
    constexpr Poly e{42, F7};
    std::cout << "e(x) = " << e << "\n";

    // ------------------------------------------------------------
    // 6. From dynamic container (truncated/padded)
    // ------------------------------------------------------------
    const std::vector<int> vec_f{9, 8, 7, 6};
    const Poly f{vec_f, F7};
    std::cout << "f(x) = " << f << "\n";

    std::cout << "\n ---------------- Operations ----------------\n";

    std::cout << "d(x) + e(x) = " << (d + e) << "\n";
    std::cout << "d(x) - e(x) = " << (d - e) << "\n";
    std::cout << "b(x) * c(x) = " << (b * c) << "\n";

    return 0;
}
