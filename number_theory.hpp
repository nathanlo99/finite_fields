
#pragma once

#include <stdexcept>
#include <string>

namespace NumberTheory {
constexpr bool is_prime_slow(int p) {
  p = std::abs(p);
  if (p < 2)
    return false;
  for (int i = 2; i * i <= p; ++i) {
    if (p % i == 0)
      return false;
  }
  return true;
}

template <class IntegerType>
constexpr IntegerType pow_mod(const IntegerType base,
                              const IntegerType exponent,
                              const IntegerType mod) {
  if (exponent < 0)
    throw std::runtime_error(
        "NumberTheory::mod_pow expects a non-negative exponent, got " +
        std::to_string(exponent));
  IntegerType result = 1;
  IntegerType pow = base;
  for (int power_of_two = 1; power_of_two < exponent; power_of_two <<= 1) {
    if (exponent & power_of_two) {
      result = (result * pow) % mod;
    }
    pow = (pow * pow) % mod;
  }
  return result;
}

template <class IntegerType>
constexpr IntegerType inv_mod(const IntegerType num, const IntegerType mod) {
  return pow_mod(num, mod - 2, mod);
}

}; // namespace NumberTheory
