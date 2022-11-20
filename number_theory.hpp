
#pragma once

#include "random.hpp"

#include <iostream>
#include <stdexcept>
#include <string>
#include <unordered_map>

namespace NumberTheory {
template <class IntegerType> constexpr bool is_prime_slow(const IntegerType p) {
  assert(p >= 0);
  if (p < 2)
    return false;
  for (IntegerType i = 2; i * i <= p; ++i) {
    if (p % i == 0)
      return false;
  }
  return true;
}

template <class IntegerType>
constexpr IntegerType add_mod(const IntegerType a, const IntegerType b,
                              const IntegerType mod) {
  assert(0 <= a && a < mod);
  assert(0 <= b && b < mod);
  return (a + b) % mod;
}

template <class IntegerType>
constexpr IntegerType sub_mod(const IntegerType a, const IntegerType b,
                              const IntegerType mod) {
  assert(0 <= a && a < mod);
  assert(0 <= b && b < mod);
  return (a + mod - b) % mod;
}

template <class IntegerType>
constexpr IntegerType mul_mod(const IntegerType a, const IntegerType b,
                              const IntegerType mod) {
  assert(0 <= a && a < mod);
  assert(0 <= b && b < mod);
  return (a * b) % mod;
}

template <class IntegerType>
constexpr IntegerType pow_mod(IntegerType base, IntegerType exponent,
                              const IntegerType mod) {
  assert(0 <= base && base < mod);
  if (exponent < 0)
    throw std::runtime_error("pow_mod() expects a non-negative exponent, got " +
                             std::to_string(exponent));
  IntegerType result = 1;
  while (exponent != 0) {
    if (exponent % 2 != 0)
      result = mul_mod(result, base, mod);
    base = mul_mod(base, base, mod);
    exponent /= 2;
  }
  return result;
}

// Returns [base ^ (2 ^ exp)] % mod
template <class IntegerType>
constexpr IntegerType pow_2_pow_mod(const IntegerType base,
                                    const IntegerType exp,
                                    const IntegerType mod) {
  IntegerType result = base;
  for (IntegerType i = 0; i < exp; ++i)
    result = mul_mod(result, result, mod);
  return result;
}

template <class IntegerType>
constexpr IntegerType inv_mod_slow(const IntegerType num,
                                   const IntegerType mod) {
  assert(0 < num && num < mod);
  return pow_mod(num, mod - 2, mod);
}

template <class IntegerType>
constexpr bool is_quadratic_residue(const IntegerType num,
                                    const IntegerType p) {
  if (!is_prime_slow(p))
    throw std::runtime_error(
        "is_quadratic_residue expects a prime modulus, got: " +
        std::to_string(p));
  if (p == 2)
    return true;
  const IntegerType fermat_value = pow_mod(num, (p - 1) / 2, p);
  return fermat_value == 0 || fermat_value == 1;
}

// Implementations of square roots mod primes taken from
// https://www.rieselprime.de/ziki/Modular_square_root
template <class IntegerType>
constexpr IntegerType sqrt_mod_8k_plus_5(const IntegerType num,
                                         const IntegerType p) {
  if (p < 0 || p % 8 != 5 || !is_prime_slow(p))
    throw std::runtime_error("sqrt_mod_8k_plus_5 expects a positive prime "
                             "modulus congruent to 5 mod 8, got " +
                             std::to_string(p));
  const IntegerType two = 2;
  const IntegerType twice_num = mul_mod(two, num, p);
  // The source above mentions (p - 5) / 8 but these are equivalent since
  // integer division rounds downwards
  const IntegerType v = pow_mod(twice_num, p / 8, p);
  const IntegerType num_v = mul_mod(num, v, p);
  const IntegerType i = mul_mod(mul_mod(twice_num, v, p), v, p);
  // i is guaranteed not to be 0, so i - 1 >= 0
  return mul_mod(num_v, i - 1, p);
}

template <class IntegerType>
IntegerType sqrt_mod_8k_plus_1(const IntegerType num, const IntegerType p) {
  if (p < 0 || p % 8 != 1 || !is_prime_slow(p))
    throw std::runtime_error("sqrt_mod_8k_plus_1 expects a positive prime "
                             "modulus congruent to 1 mod 8, got " +
                             std::to_string(p));

  IntegerType e = 0, q = p - 1;
  while (q % 2 == 0) {
    q /= 2;
    e += 1;
  }

  IntegerType z = 0;
  do {
    const IntegerType x = util::random_int64_t(2, p - 1);
    z = pow_mod(x, q, p);
  } while (pow_2_pow_mod(z, e - 1, p) == 1);

  IntegerType y = z, r = e, x = pow_mod(num, (q - 1) / 2, p),
              v = mul_mod(num, x, p), w = mul_mod(v, x, p);
  while (w != 1) {
    IntegerType pow = w;
    IntegerType k = 0;
    while (pow != 1) {
      pow = mul_mod(pow, pow, p);
      k += 1;
    }
    const IntegerType d = pow_2_pow_mod(y, r - k - 1, p);
    y = mul_mod(d, d, p);
    r = k;
    v = mul_mod(d, v, p);
    w = mul_mod(w, y, p);
  }
  return v;
}

template <class IntegerType>
IntegerType sqrt_mod(const IntegerType num, const IntegerType p) {
  if (p < 0 || !is_prime_slow(p))
    throw std::runtime_error("sqrt_mod expects a positive prime modulus, got " +
                             std::to_string(p));
  if (num < 0 || num >= p)
    throw std::runtime_error(
        "sqrt_mod expects 0 <= num < p, got num = " + std::to_string(num) +
        ", p = " + std::to_string(p));

  if (p == 2 || num == 0)
    return num;
  else if (p % 4 == 3)
    return pow_mod(num, (p + 1) / 4, p);
  else if (p % 8 == 5) {
    return sqrt_mod_8k_plus_5(num, p);
  } else {
    return sqrt_mod_8k_plus_1(num, p);
  }
}

}; // namespace NumberTheory
