
#pragma once

#include "field.hpp"
#include "fraction.hpp"

#include <cstdint>
#include <stdexcept>

template <class Fraction = Fraction<>>
class RationalField : AbstractField<Fraction> {
public:
  // NOTE: These are the two types which have to be exposed as a subclass of
  // AbstractField
  using value_t = Fraction;
  using element_t = FieldElement<RationalField>;

public:
  constexpr inline element_t operator()(const int64_t num) const {
    return element_t(value_t(num, 1), *this);
  }
  constexpr inline int64_t characteristic() const override { return 0; }

  constexpr inline value_t zero() const override { return Fraction(0, 1); }
  constexpr inline value_t one() const override { return Fraction(1, 1); }

  constexpr inline value_t neg(const value_t a) const override { return -a; }
  constexpr inline value_t add(const value_t a,
                               const value_t b) const override {
    return a + b;
  }
  constexpr inline value_t sub(const value_t a,
                               const value_t b) const override {
    return a - b;
  }
  constexpr inline value_t inv(const value_t a) const override {
    return a.inv();
  }
  constexpr inline value_t mul(const value_t a,
                               const value_t b) const override {
    return a * b;
  }
  constexpr inline value_t div(const value_t a,
                               const value_t b) const override {
    return a / b;
  }
  constexpr inline bool eq(const value_t a, const value_t b) const override {
    return a == b;
  }

  friend std::ostream &operator<<(std::ostream &os,
                                  const RationalField &field) {
    return os << "Rational field";
  }
};