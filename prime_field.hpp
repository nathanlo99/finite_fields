
#pragma once

#include "error.hpp"
#include "field.hpp"
#include "number_theory.hpp"

#include <cassert>
#include <cmath>
#include <cstdint>
#include <stdexcept>

namespace nt = NumberTheory;

template <class IntegerType = int64_t>
class PrimeField : AbstractField<IntegerType> {
public:
  using value_t = IntegerType;
  using element_t = FieldElement<PrimeField>;

  const value_t p;

  constexpr PrimeField(const value_t p) : p(p) {
    if (p < 2)
      throw math_error("PrimeField modulus must be at least 2");
    if (std::sqrt(std::numeric_limits<IntegerType>::max()) < p)
      throw math_error()
          << "PrimeField integer type was not large enough to support modulus "
          << p;
    if (!nt::is_prime(p))
      throw math_error() << "Modulus for PrimeField must be prime: " << p;
  }

  constexpr inline void assert_in_bounds(const value_t num) const {
    assert(0 <= num && num < p);
  }

  constexpr inline element_t operator()(const int64_t num) const {
    return element_t(integer(num), *this);
  }

  constexpr inline int64_t characteristic() const override { return p; }

  constexpr inline value_t zero() const override { return 0; }
  constexpr inline value_t one() const override { return 1; }

  constexpr inline value_t integer(const int64_t value) const override {
    return ((value % p) + p) % p;
  }
  constexpr inline element_t element(const value_t value) const {
    return element_t(integer(value), *this);
  }

  constexpr inline value_t neg(const value_t a) const override {
    assert_in_bounds(a);
    return a == 0 ? 0 : p - a;
  }
  constexpr inline value_t add(const value_t a,
                               const value_t b) const override {
    assert_in_bounds(a);
    assert_in_bounds(b);
    return nt::add_mod<value_t>(a, b, p);
  }
  constexpr inline value_t sub(const value_t a,
                               const value_t b) const override {
    assert_in_bounds(a);
    assert_in_bounds(b);
    return nt::sub_mod<value_t>(a, b, p);
  }

  constexpr inline value_t inv(const value_t a) const override {
    assert_in_bounds(a);
    if (a == 0)
      throw math_error("Division by zero");
    return nt::inv_mod<value_t>(a, p);
  }
  constexpr inline value_t mul(const value_t a,
                               const value_t b) const override {
    assert_in_bounds(a);
    assert_in_bounds(b);
    return nt::mul_mod<value_t>(a, b, p);
  }
  constexpr inline value_t div(const value_t a,
                               const value_t b) const override {
    assert_in_bounds(a);
    assert_in_bounds(b);
    return mul(a, inv(b));
  }

  constexpr inline value_t pow(const value_t a,
                               const uint64_t exp) const override {
    assert_in_bounds(a);
    return nt::pow_mod(a, exp, p);
  }

  constexpr inline bool eq(const value_t a, const value_t b) const override {
    assert_in_bounds(a);
    assert_in_bounds(b);
    return a == b;
  }

  constexpr inline bool operator==(const PrimeField &other) const {
    return p == other.p;
  }

public:
  friend std::ostream &operator<<(std::ostream &os, const PrimeField &field) {
    return os << "Finite Field of size " << field.p;
  }
};
