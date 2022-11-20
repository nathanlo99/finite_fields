
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

// Uses the extended Euclidean algorithm to find integer solutions to
// num * inv + k * p = 1
template <class IntegerType>
constexpr IntegerType inv_mod(const IntegerType num, const IntegerType p) {
  assert(0 < num && num < p);
  IntegerType r0 = p, r1 = num, s0 = 1, s1 = 0, t0 = 0, t1 = 1;
  while (r1 != 0) {
    const IntegerType q = r0 / r1, r2 = r0 % r1;
    std::tie(r0, s0, t0, r1, s1, t1) =
        std::make_tuple(r1, s1, t1, r2, s0 - s1 * q, t0 - t1 * q);
  }
  return t0 >= 0 ? t0 : t0 + p;
}

template <class IntegerType>
constexpr IntegerType jacobi_symbol(IntegerType a, IntegerType n) {
  if (n <= 0 || n % 2 == 0)
    throw std::runtime_error(
        "jacobi_symbol expected n to be positive and odd, but got " +
        std::to_string(n));
  IntegerType result = 1;
  while (true) {
    a %= n;
    if (a == 0)
      return (n == 1) ? result : 0;
    IntegerType h = 0;
    while (a % 2 == 0) {
      a /= 2;
      h++;
    }
    const IntegerType remainder = n % 8;
    if (h % 2 != 0 && remainder != 1 && remainder != 7)
      result = -result;
    if (a % 4 != 1 && remainder != 1 && remainder != 5)
      result = -result;
    std::swap(a, n);
  }
}

template <class IntegerType>
constexpr bool is_quadratic_residue(const IntegerType num,
                                    const IntegerType p) {
  assert(0 <= num && num < p);
  if (!is_prime_slow(p))
    throw std::runtime_error(
        "is_quadratic_residue expects a prime modulus, got: " +
        std::to_string(p));
  if (p == 2 || num == 0)
    return true;
  return jacobi_symbol(num, p) == 1;
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
