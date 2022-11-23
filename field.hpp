
#pragma once

#include "error.hpp"

#include <iostream>
#include <stdexcept>

template <class Value> class AbstractField {
public:
  using value_t = Value;

  virtual int64_t characteristic() const = 0;

  virtual value_t zero() const = 0;
  virtual value_t one() const = 0;

  // The default is a naive generic double-and-add, but subclasses can override
  // this with custom behaviour, see PrimeField::integer()
  virtual value_t integer(int64_t value) const {
    bool negative = false;
    if (value < 0) {
      value = -value;
      negative = true;
    }

    value_t result = zero(), pow = one();
    for (int power_of_two = 1; power_of_two < value; power_of_two <<= 1) {
      if (value & power_of_two)
        result = add(result, pow);
      pow = add(pow, pow);
    }
    return negative ? neg(result) : result;
  }

  virtual value_t neg(const value_t a) const = 0;
  virtual value_t add(const value_t a, const value_t b) const = 0;
  virtual value_t sub(const value_t a, const value_t b) const {
    return add(a, neg(b));
  }

  virtual value_t inv(const value_t a) const = 0;
  virtual value_t mul(const value_t a, const value_t b) const = 0;
  virtual value_t div(const value_t a, const value_t b) const {
    return mul(a, inv(b));
  }
  virtual value_t pow(const value_t a, const uint64_t exp) const {
    value_t base = a, result = one();
    uint64_t expt = exp;
    while (expt > 0) {
      if (expt % 2 == 1)
        result = mul(result, base);
      base = mul(base, base);
    }
    return result;
  }

  virtual bool eq(const value_t a, const value_t b) const = 0;

  bool operator==(const AbstractField &other) { return false; }
};

template <class Field> struct FieldElement {
  using element_t = FieldElement<Field>;
  using value_t = typename Field::value_t;

  const Field &field;
  value_t value;

  FieldElement(const value_t value, const Field &field)
      : field(field), value(value) {}
  FieldElement(const FieldElement &other)
      : field(other.field), value(other.value) {}
  FieldElement &operator=(const FieldElement &other) {
    if (field != other.field)
      throw math_error("Field element cannot be assigned to field "
                       "element with different base field");
    value = other.value;
    return *this;
  }

  constexpr element_t zero(const Field &field) const {
    return FieldElement(field.zero(), field);
  }
  constexpr element_t one(const Field &field) const {
    return FieldElement(field.one(), field);
  }

  constexpr element_t operator-() const {
    return element_t(field.neg(value), field);
  }
  constexpr element_t operator+(const element_t &other) const {
    return element_t(field.add(value, other.value), field);
  }
  constexpr element_t operator-(const element_t &other) const {
    return element_t(field.sub(value, other.value), field);
  }
  constexpr element_t &operator+=(const element_t &other) {
    return *this = *this + other;
  }
  constexpr element_t &operator-=(const element_t &other) {
    return *this = *this - other;
  }

  constexpr element_t inv() const { return element_t(field.inv(value), field); }
  constexpr element_t operator*(const element_t &other) const {
    return element_t(field.mul(value, other.value), field);
  }
  constexpr element_t operator/(const element_t &other) const {
    return element_t(field.div(value, other.value), field);
  }
  constexpr element_t operator^(const uint64_t exp) const {
    return element_t(field.pow(value, exp), field);
  }

  constexpr bool operator==(const element_t &other) const {
    return field.eq(value, other.value);
  }

  friend element_t operator*(const int a, const element_t &b) {
    return b.field(a) * b;
  }
  friend element_t operator*(const element_t &a, const int b) {
    return a * a.field(b);
  }

  friend bool operator==(const int a, const element_t &b) {
    return b.field.eq(a, b.value);
  }
  friend bool operator==(const element_t &a, const int b) {
    return a.field.eq(a.value, b);
  }

  friend std::ostream &operator<<(std::ostream &os, const FieldElement &val) {
    return os << val.value;
  }
};
