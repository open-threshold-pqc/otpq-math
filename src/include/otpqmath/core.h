#pragma once

/**
 * @file core.h
 * @brief Unified include point for all algebraic structures in OTPQ-Math.
 *
 * This header aggregates the most commonly used algebra backends and
 * polynomial types so that users can access the full algebraic subsystem
 * via a single include.
 *
 * OTPQ-Math separates *concepts* (abstract algebraic interfaces) from
 * *implementations* (modular integer rings, finite fields, complex fields,
 * polynomial rings, etc.).  This header pulls together the concrete
 * implementations that model the concepts defined in:
 *
 *      otpq::math::als::core
 *
 * It includes:
 *
 * ────────────────────────────────────────────────────────────────
 *  RINGS
 * ────────────────────────────────────────────────────────────────
 *   • ModRing
 *        Ring of integers modulo n (ℤ/nℤ), supporting add/sub/mul/neg.
 *
 * ────────────────────────────────────────────────────────────────
 *  FIELDS
 * ────────────────────────────────────────────────────────────────
 *   • PrimeField
 *        Finite field 𝔽ₚ, with full inverse/division support.
 *
 *   • ComplexField
 *        Floating-point complex numbers ℂ, used in FFT/NTT and lattice
 *        sampling (e.g., Falcon, GPV).
 *
 * ────────────────────────────────────────────────────────────────
 *  POLYNOMIALS
 * ────────────────────────────────────────────────────────────────
 *   • Polynomial
 *        Fixed-size polynomial container (degree < N) with efficient
 *        data access (operator[], raw(), data()).
 *
 *   • PolynomialRing
 *        Quotient ring  Alg[x] / M(x), parameterized by:
 *            – coefficient algebra (e.g., PrimeField)
 *            – polynomial modulus (e.g., xⁿ−1, xⁿ+1)
 *
 * All structures are:
 *   - constexpr-friendly,
 *   - lightweight (no dynamic allocations unless explicitly requested),
 *   - usable in SIMD, FFT/NTT, or hardware-accelerated backends,
 *   - compatible with compile-time algebraic concepts (Ring, Field).
 *
 * @note This header is intentionally broad. For performance-critical code
 *       you may prefer including only the specific algebra modules needed.
 */

#include <otpqmath/als/rings/numbers/mod_ring.h>
#include <otpqmath/als/fields/numbers/prime_field.h>
#include <otpqmath/als/fields/numbers/complex_field.h>
#include <otpqmath/als/rings/poly/polynomial.h>
#include <otpqmath/als/rings/poly/polynomial_ring.h>
