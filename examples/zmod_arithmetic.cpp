#include <iostream>
#include <zmod_arithmetic.h>

namespace mod = otpq::math::als::arith;

int main() {
    std::cout << "=== ZmodNRing Example (ℤ/nℤ) ===\n";

    const mod::ZmodNRing R7{7};

    std::cout << "modulus = " << R7.mod << "\n";
    std::cout << "3 + 6 = " << R7.add(3, 6) << "   (mod 7)\n";
    std::cout << "4 * 5 = " << R7.mul(4, 5) << "   (mod 7)\n";
    std::cout << "2 - 5 = " << R7.sub(2, 5) << "   (mod 7)\n";
    std::cout << "-3    = " << R7.neg(3) << "     (mod 7)\n";

    std::cout << "\n=== ZmodNField Example (𝔽₇) ===\n";


    const mod::ZmodNField F7{7};

    std::cout << "modulus = " << F7.mod << "\n";
    std::cout << "3 + 6  = " << F7.add(3, 6) << "\n";
    std::cout << "4 * 5  = " << F7.mul(4, 5) << "\n";
    std::cout << "inv(3) = " << F7.inv(3) << "\n";
    std::cout << "2 / 4  = " << F7.div(2, 4) << "\n";

    std::cout << "\nMultiplicative inverses in 𝔽₇:\n";
    for (int i = 1; i < 7; ++i)
        std::cout << "  inv(" << i << ") = " << F7.inv(i) << "\n";

    return 0;
}
