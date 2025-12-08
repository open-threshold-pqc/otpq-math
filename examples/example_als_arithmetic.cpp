#include <iostream>
#include <complex>

#include <otpqmath/als/ring/numbers/mod_ring.h>
#include <otpqmath/als/field/numbers/prime_field.h>
#include <otpqmath/als/field/numbers/complex_field.h>

using otpq::math::als::ring::numbers::ModRing;
using otpq::math::als::field::number::PrimeField;
using otpq::math::als::field::number::ComplexField;

int main() {
    std::cout << "=== ModRing Example  (ℤ/nℤ) ===\n";

    constexpr ModRing<int> R7{7};

    std::cout << "modulus = " << R7.mod << "\n";
    std::cout << "3 + 6 = " << R7.add(3, 6) << "   (mod 7)\n";
    std::cout << "4 * 5 = " << R7.mul(4, 5) << "   (mod 7)\n";
    std::cout << "2 - 5 = " << R7.sub(2, 5) << "   (mod 7)\n";
    std::cout << "-3    = " << R7.neg(3) << "     (mod 7)\n";


    // ----------------------------------------------------------------------
    //  Finite Field 𝔽₇
    // ----------------------------------------------------------------------
    std::cout << "\n=== PrimeField Example (𝔽₇) ===\n";

    constexpr PrimeField<int> F7{7};

    std::cout << "modulus = " << F7.mod << "\n";
    std::cout << "3 + 6  = "  << F7.add(3, 6) << "\n";
    std::cout << "4 * 5  = "  << F7.mul(4, 5) << "\n";
    std::cout << "inv(3) = "  << F7.inv(3) << "\n";
    std::cout << "2 / 4  = "  << F7.div(2, 4) << "\n";

    std::cout << "\nMultiplicative inverses in 𝔽₇:\n";
    for (int i = 1; i < 7; ++i)
        std::cout << "  inv(" << i << ") = " << F7.inv(i) << "\n";


    // ----------------------------------------------------------------------
    //  Complex Field ℂ
    // ----------------------------------------------------------------------
    std::cout << "\n=== ComplexField Example (ℂ) ===\n";

    ComplexField<double> C;

    std::complex<double> a{1.0, 2.0};
    std::complex<double> b{3.0, -1.0};

    std::cout << "a        = " << a << "\n";
    std::cout << "b        = " << b << "\n";
    std::cout << "a + b    = " << C.add(a, b) << "\n";
    std::cout << "a - b    = " << C.sub(a, b) << "\n";
    std::cout << "a * b    = " << C.mul(a, b) << "\n";
    std::cout << "-a       = " << C.neg(a) << "\n";
    std::cout << "inv(a)   = " << C.inv(a) << "\n";
    std::cout << "a / b    = " << C.div(a, b) << "\n";

    return 0;
}
