
#pragma once

#include "field.hpp"
#include "number_theory.hpp"

#include <cmath>
#include <cstdint>
#include <stdexcept>

namespace nt = NumberTheory;

template <class IntegerType = int64_t>
class PrimeField : AbstractField<IntegerType> {
public:
  using value_t = IntegerType;
  using element_t = FieldElement<PrimeField<IntegerType>>;

  const value_t p;

  constexpr PrimeField(const value_t p) : p(p) {
    if (std::sqrt(std::numeric_limits<IntegerType>::max()) < p)
      throw std::runtime_error(
          "PrimeField integer type was not large enough to support modulus " +
          std::to_string(p));
    if (!nt::is_prime_slow(p))
      throw std::runtime_error("Modulus for PrimeField must be prime: " +
                               std::to_string(p));
  }

  constexpr inline element_t value(const value_t num) const {
    return element_t(val(num), this);
  }
  constexpr inline value_t zero() const override { return 0; }
  constexpr inline value_t one() const override { return 1; }
  constexpr inline value_t integer(const int value) const override {
    return val(value);
  }

  constexpr inline value_t neg(const value_t a) const override {
    return val(p - a);
  }
  constexpr inline value_t add(const value_t a,
                               const value_t b) const override {
    return val(a + b);
  }
  constexpr inline value_t sub(const value_t a,
                               const value_t b) const override {
    return val(a - b);
  }

  constexpr inline value_t inv(const value_t a) const override {
    if (a == 0)
      throw std::runtime_error("Division by zero");
    return nt::inv_mod<value_t>(a, p);
  }
  constexpr inline value_t mul(const value_t a,
                               const value_t b) const override {
    return val(a * b);
  }
  constexpr inline value_t div(const value_t a,
                               const value_t b) const override {
    if (b == 0)
      throw std::runtime_error("Division by zero");
    return val(a * inv(b));
  }

  constexpr inline bool eq(const value_t a, const value_t b) const override {
    return a == b;
  }

private:
  constexpr inline value_t val(const value_t value) const {
    return static_cast<value_t>(((value % p) + p) % p);
  }

public:
  friend std::ostream &operator<<(std::ostream &os, const PrimeField &field) {
    return os << "Finite Field of size " << field.p;
  }
};
