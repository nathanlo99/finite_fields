
#pragma once

#include "random.hpp"

#include <iostream>
#include <stdexcept>
#include <string>

namespace NumberTheory {
constexpr bool is_prime_slow(int p) {
  assert(p >= 0);
  if (p < 2)
    return false;
  for (int i = 2; i * i <= p; ++i) {
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
  return (((a - b) % mod) + mod) % mod;
}

template <class IntegerType>
constexpr IntegerType mul_mod(const IntegerType a, const IntegerType b,
                              const IntegerType mod) {
  assert(0 <= a && a < mod);
  assert(0 <= b && b < mod);
  return (a * b) % mod;
}

template <class IntegerType>
constexpr IntegerType pow_mod(const IntegerType base,
                              const IntegerType exponent,
                              const IntegerType mod) {
  assert(0 <= base && base < mod);
  if (exponent < 0)
    throw std::runtime_error("pow_mod() expects a non-negative exponent, got " +
                             std::to_string(exponent));
  IntegerType result = 1;
  IntegerType pow = base;
  for (int power_of_two = 1; power_of_two <= exponent; power_of_two <<= 1) {
    if (exponent & power_of_two) {
      result = mul_mod(result, pow, mod);
    }
    pow = mul_mod(pow, pow, mod);
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
constexpr IntegerType inv_mod(const IntegerType num, const IntegerType mod) {
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
  return p == 2 || pow_mod(num, (p - 1) / 2, p) == 1;
}

// Implementation from https://www.rieselprime.de/ziki/Modular_square_root
template <class IntegerType>
constexpr IntegerType sqrt_mod(IntegerType num, IntegerType p) {
  if (!is_prime_slow(p))
    throw std::runtime_error("sqrt_mod expects a prime modulus, got " +
                             std::to_string(p));
  num %= p;
  if (p == 2 || num == 0)
    return num;
  else if (p % 4 == 3)
    return pow_mod(num, (p + 1) / 4, p);
  else if (p % 8 == 5) {
    const IntegerType two = 2;
    const IntegerType twice_num = mul_mod(two, num, p);
    // The source above mentions (p - 5) / 8 but these are equivalent since
    // integer division rounds downwards
    const IntegerType v = pow_mod(twice_num, p / 8, p);
    const IntegerType num_v = mul_mod(num, v, p);
    const IntegerType i = mul_mod(mul_mod(twice_num, v, p), v, p);
    // i is guaranteed not to be 0, so i - 1 >= 0
    return mul_mod(num_v, i - 1, p);
  } else {
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
}

}; // namespace NumberTheory
