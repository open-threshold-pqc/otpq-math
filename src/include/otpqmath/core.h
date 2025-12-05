#pragma once

/**
 * @file core.h
 * @brief Unified include point for algebraic structures (rings and fields).
 *
 * This module provides the algebra backends used throughout OTPQ-Math for:
 *   - polynomial arithmetic,
 *   - lattice-based cryptography,
 *   - number theory,
 *   - FFT/NTT computations.
 *
 * The types included here implement the algebraic concepts defined in
 * `otpq::math::als::core` (Ring, Field), enabling generic and efficient
 * coefficient arithmetic for polynomial types.
 *
 * Included categories:
 *
 *  ─ Rings:
 *      • ModRing
 *          Integer residue ring ℤ/nℤ.
 *
 *  ─ Fields:
 *      • PrimeField
 *          Finite field 𝔽ₚ for prime modulus p.
 *
 *      • ComplexField
 *          Field of complex numbers ℂ, used for FFT/NTT and
 *          floating-point lattice operations (e.g., Falcon).
 *
 * All structures are lightweight, constexpr-friendly, and suitable for
 * SIMD or hardware-accelerated arithmetic via backend dispatch.
 */

#include <otpqmath/numbers/rings/mod_ring.h>
#include <otpqmath/numbers/fields/prime_field.h>
#include <otpqmath/numbers/fields/complex_field.h>
