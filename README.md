<div align="center">

<img src="assets/logo.png" alt="OTPQ-MATH Logo" width="180" />

# **OTPQ-MATH**
### High-performance algebra, polynomial fields, and modular arithmetic for C++

---

<!-- Badges -->
<p>
  <img src="https://img.shields.io/github/languages/top/yourname/otpq-math" />
  <img src="https://img.shields.io/github/last-commit/yourname/otpq-math" />
  <img src="https://img.shields.io/github/license/yourname/otpq-math" />
</p>

<h4>
  <a href="#eyes-usage">Usage</a>
  <span> • </span>
  <a href="#toolbox-getting-started">Getting Started</a>
  <span> • </span>
  <a href="#compass-roadmap">Roadmap</a>
  <span> • </span>
  <a href="#wave-contributing">Contribute</a>
</h4>

</div>

---

# :star2: About OTPQ-MATH

**OTPQ-MATH** is a modern C++ library implementing:

- Modular arithmetic (Z mod q) over rings and fields
- Fixed-degree polynomial arithmetic
- Polynomial fields / quotient rings
- Cyclotomic and anti-cyclotomic moduli
- Support for lattice-cryptography-style polynomial operations
- Extensible algebraic structures (ALS)

It is designed to be **fast, small, header-only**, and suitable for:

- Cryptography (RLWE, NTT, Ring-LWE, etc.)
- Computational algebra
- Numeric research
- Polynomial field experimentation

---

# :space_invader: Features

- ⚙️ **ZmodN ring/field algebra** (`Als` abstraction)
- 🧮 **Polynomial<N>** fixed-size polynomial class
- 🔁 **PolynomialField** supporting:
    - Addition, subtraction, multiplication
    - Modular reduction
    - **Multiplicative inverse** via Extended Euclidean Algorithm
- 🌀 **Cyclotomic modulus (xⁿ − 1)**
- 🔄 **Anti-cyclotomic modulus (xⁿ + 1)**
- ❌ Safe rejection of inversion for `NoPolyModulus`
- 🧱 Extensible design for arbitrary **M(x)** in the future
- ⚡ Write-once header-only architecture
- 🧪 Unit-test friendly structure

---

# :toolbox: Getting Started

## :bangbang: Requirements

- **C++20** or later
- A modern compiler:
    - GCC ≥ 11
    - Clang ≥ 13
    - MSVC ≥ 2022
- CMake ≥ 3.20 (optional)

---

# :gear: Installation

### Option 1 — Header-only (recommended)

Just include the library directory:

```cpp
#include "polynomial_fields.h"
#include "zmod_arithmetic.h"
#include "polynomial.h"
