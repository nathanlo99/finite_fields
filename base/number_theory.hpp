
#pragma once

#include "error.hpp"
#include "random.hpp"

#include <array>
#include <iostream>
#include <stdexcept>
#include <string>
#include <unordered_map>

namespace NumberTheory {

template <class IntegerType>
constexpr inline IntegerType gcd(IntegerType a, IntegerType b) {
  while (b != 0)
    std::tie(a, b) = std::make_pair(b, a % b);
  return a;
}

// Returns a pair (s, k) such that num = 2^s * k, s >= 0, and k is odd
template <class IntegerType>
constexpr std::pair<IntegerType, IntegerType>
odd_factorization(IntegerType num) {
  IntegerType s = 0;
  while (num % 2 == 0) {
    num /= 2;
    s++;
  }
  return std::make_pair(s, num);
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
    throw math_error() << "pow_mod() expects a non-negative exponent, got "
                       << exponent;
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
    throw math_error()
        << "jacobi_symbol expected n to be positive and odd, but got " << n;
  IntegerType result = 1;
  while (true) {
    a %= n;
    if (a == 0)
      return (n == 1) ? result : 0;
    const auto &[h, a1] = odd_factorization(a);
    const IntegerType remainder = n % 8;
    if (h % 2 != 0 && remainder != 1 && remainder != 7)
      result = -result;
    if (a1 % 4 != 1 && remainder != 1 && remainder != 5)
      result = -result;
    a = n;
    n = a1;
  }
}

// Given n, s, d, a with n - 1 = 2^s * d, returns true if and only if a is a
// Miller-Rabin witness for n
template <class IntegerType>
constexpr bool is_miller_rabin_witness(const IntegerType n, IntegerType s,
                                       const IntegerType d,
                                       const IntegerType a) {
  assert(odd_factorization(n - 1) == std::make_pair(s, d));
  IntegerType x = pow_mod(a, d, n), y = 0;
  while (s > 0) {
    y = mul_mod(x, x, n);
    if (y == 1 && x != 1 && x != n - 1)
      return false;
    x = y;
    --s;
  }
  return y == 1;
}

template <class IntegerType>
constexpr bool is_likely_prime_miller_rabin(const IntegerType n,
                                            const IntegerType s,
                                            const IntegerType d,
                                            const IntegerType num_iterations) {
  for (IntegerType iter = 0; iter < num_iterations; ++iter) {
    const IntegerType a = util::random_int64_t(2, n - 1);
    if (!is_miller_rabin_witness(n, s, d, a))
      return false;
  }
  return true;
}

template <class IntegerType>
constexpr bool is_prime_miller_rabin(const IntegerType n) {
  if (n == 2 || n == 3)
    return true;
  if (n < 2 || n % 2 == 0 || n % 3 == 0)
    return false;
  const auto &[s, d] = odd_factorization(n - 1);
  if (n < 1373653LL)
    return is_miller_rabin_witness<IntegerType>(n, s, d, 2) &&
           is_miller_rabin_witness<IntegerType>(n, s, d, 3);
  if (n < 9080191LL)
    return is_miller_rabin_witness<IntegerType>(n, s, d, 31) &&
           is_miller_rabin_witness<IntegerType>(n, s, d, 73);
  if (n < 4759123141LL)
    return is_miller_rabin_witness<IntegerType>(n, s, d, 2) &&
           is_miller_rabin_witness<IntegerType>(n, s, d, 7) &&
           is_miller_rabin_witness<IntegerType>(n, s, d, 61);
  if (n < 1122004669633LL)
    return is_miller_rabin_witness<IntegerType>(n, s, d, 2) &&
           is_miller_rabin_witness<IntegerType>(n, s, d, 13) &&
           is_miller_rabin_witness<IntegerType>(n, s, d, 23) &&
           is_miller_rabin_witness<IntegerType>(n, s, d, 1662803);
  if (n < 2152302898747LL)
    return is_miller_rabin_witness<IntegerType>(n, s, d, 2) &&
           is_miller_rabin_witness<IntegerType>(n, s, d, 3) &&
           is_miller_rabin_witness<IntegerType>(n, s, d, 5) &&
           is_miller_rabin_witness<IntegerType>(n, s, d, 7) &&
           is_miller_rabin_witness<IntegerType>(n, s, d, 11);
  if (n < 3474749660383LL)
    return is_miller_rabin_witness<IntegerType>(n, s, d, 2) &&
           is_miller_rabin_witness<IntegerType>(n, s, d, 3) &&
           is_miller_rabin_witness<IntegerType>(n, s, d, 5) &&
           is_miller_rabin_witness<IntegerType>(n, s, d, 7) &&
           is_miller_rabin_witness<IntegerType>(n, s, d, 11) &&
           is_miller_rabin_witness<IntegerType>(n, s, d, 13);
  if (n < 341550071728321LL)
    return is_miller_rabin_witness<IntegerType>(n, s, d, 2) &&
           is_miller_rabin_witness<IntegerType>(n, s, d, 3) &&
           is_miller_rabin_witness<IntegerType>(n, s, d, 5) &&
           is_miller_rabin_witness<IntegerType>(n, s, d, 7) &&
           is_miller_rabin_witness<IntegerType>(n, s, d, 11) &&
           is_miller_rabin_witness<IntegerType>(n, s, d, 13) &&
           is_miller_rabin_witness<IntegerType>(n, s, d, 17);
  return is_likely_prime_miller_rabin<IntegerType>(n, s, d, 100);
}

template <class IntegerType>
constexpr inline bool is_prime(const IntegerType n) {
  assert(n >= 0);
  return is_prime_miller_rabin(n);
}

template <class IntegerType>
constexpr inline IntegerType next_prime(IntegerType n) {
  if (n < 0)
    throw math_error() << "next_prime() expected a non-negative integer, got "
                       << n;

  // The numbers coprime to 30 are 1 7 11 13 17 19 23 29
  constexpr std::array<IntegerType, 8> small_primes = {2, 2, 3, 5, 5, 7, 7, 11};
  constexpr std::array<IntegerType, 30> offset_table = {
      1, 6, 5, 4, 3, 2, //  0 -  5
      1, 4, 3, 2, 1, 2, //  6 - 11
      1, 4, 3, 2, 1, 2, // 12 - 17
      1, 4, 3, 2, 1, 6, // 18 - 23
      5, 4, 3, 2, 1, 2, // 24 - 29
  };
  if (n < 8) [[unlikely]]
    return small_primes[n];
  while (true) {
    n += offset_table[n % 30];
    if (is_prime(n))
      return n;
  }
}

template <class IntegerType>
constexpr bool is_quadratic_residue(const IntegerType num,
                                    const IntegerType p) {
  assert(0 <= num && num < p);
  if (!is_prime(p))
    throw math_error() << "is_quadratic_residue expects a prime modulus, got: "
                       << p;
  if (p == 2 || num == 0)
    return true;
  return jacobi_symbol(num, p) == 1;
}

// Implementations of square roots mod primes taken from
// https://www.rieselprime.de/ziki/Modular_square_root
template <class IntegerType>
constexpr IntegerType sqrt_mod_8k_plus_5(const IntegerType num,
                                         const IntegerType p) {
  if (p < 0 || p % 8 != 5 || !is_prime(p))
    throw math_error() << "sqrt_mod_8k_plus_5 expects a positive prime "
                          "modulus congruent to 5 mod 8, got "
                       << p;
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
  if (p < 0 || p % 8 != 1 || !is_prime(p))
    throw math_error() << "sqrt_mod_8k_plus_1 expects a positive prime modulus "
                          "congruent to 1 mod 8, got "
                       << p;

  const auto &[e, q] = odd_factorization(p - 1);
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
  if (p < 0 || !is_prime(p))
    throw math_error() << "sqrt_mod expects a positive prime modulus, got "
                       << p;
  if (num < 0 || num >= p)
    throw math_error() << "sqrt_mod expects 0 <= num < p, got num = " << num
                       << ", p = " << p;

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
