#include <iostream>
#include <array>

#include <polynomial_fields.h>
#include <ntt.h>

using namespace otpq::math::als::polynomial;
using namespace otpq::math::als::polynomial::transforms::ntt;


int main() {

    // Choose parameters
    constexpr std::size_t N = 256;
    int modulus = 8380417;

    PolynomialField<N> p({1, 2, 3, 4, 3, 2, 1, 51}, modulus);
    PolynomialField<N> q({41, 2, 3, 4, 3, 2, 1, 51}, modulus);

    std::cout << "p(x): " << p << std::endl;
    std::cout << "q(x): " << q << std::endl;


    int primitive_root = 1753;

    auto p_ntt = toNTT(p, primitive_root);
    auto q_ntt = toNTT(q, primitive_root);

    std::cout << p_ntt << '\n';
    std::cout << q_ntt << '\n';

    p = fromNTT(p_ntt, 731434);

    std::cout << (mulNTT(p_ntt, q_ntt))  << '\n';

    return 0;
}
