
#pragma once

#include "error.hpp"
#include "number_theory.hpp"

#include <cstdint>
#include <stdexcept>

namespace nt = NumberTheory;

template <class IntegerType = int64_t> struct Fraction {
  using value_t = IntegerType;

  value_t num, denom;

  constexpr Fraction(value_t num, value_t denom = 1) {
    if (denom == 0)
      throw math_error("Fraction: division by zero");
    const bool negative = (num < 0) ^ (denom < 0);
    num = std::abs(num);
    denom = std::abs(denom);
    const value_t gcd = nt::gcd(num, denom);
    this->num = negative ? -num / gcd : num / gcd;
    this->denom = denom / gcd;
  }

  constexpr inline Fraction operator-() const { return Fraction(-num, denom); }

  friend constexpr inline Fraction operator+(const Fraction &a,
                                             const Fraction &b) {
    return Fraction(a.num * b.denom + a.denom * b.num, a.denom * b.denom);
  }
  friend constexpr inline Fraction operator-(const Fraction &a,
                                             const Fraction &b) {
    return a + (-b);
  }

  constexpr inline Fraction inv() const {
    if (num == 0)
      throw math_error("Fraction::inv(): division by zero");
    return Fraction(denom, num);
  }
  friend constexpr inline Fraction operator*(const Fraction &a,
                                             const Fraction &b) {
    return Fraction(a.num * b.num, a.denom * b.denom);
  }
  friend constexpr inline Fraction operator/(const Fraction &a,
                                             const Fraction &b) {
    return Fraction(a.num * b.denom, a.denom * b.num);
  }

  friend constexpr inline bool operator==(const Fraction &a,
                                          const Fraction &b) {
    return a.num == b.num && a.denom == b.denom;
  }

  friend std::ostream &operator<<(std::ostream &os, const Fraction &fraction) {
    if (fraction.denom == 1)
      return os << fraction.num;
    return os << fraction.num << "/" << fraction.denom;
  }
};
