
#pragma once

#include <iostream>

template <class Value> class AbstractField {
public:
  using value_t = Value;

  virtual value_t zero() const = 0;
  virtual value_t one() const = 0;

  // The default is a naive generic double-and-add, but subclasses can override
  // this with custom behaviour, see PrimeField::integer()
  virtual value_t integer(const int64_t value) const {
    value_t result = zero(), pow = one();
    for (int power_of_two = 1; power_of_two < value; power_of_two <<= 1) {
      if (value & power_of_two)
        result = add(result, pow);
      pow = add(pow, pow);
    }
    return result;
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

  virtual bool eq(const value_t a, const value_t b) const = 0;
};

template <class Field> struct FieldElement {
  using element_t = FieldElement<Field>;
  using value_t = typename Field::value_t;

  value_t value;
  const Field &field;

  FieldElement(const value_t value, const Field &field)
      : value(value), field(field) {}

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

  constexpr element_t inv() const { return element_t(field.inv(value), field); }
  constexpr element_t operator*(const element_t &other) const {
    return element_t(field.mul(value, other.value), field);
  }
  constexpr element_t operator/(const element_t &other) const {
    return element_t(field.div(value, other.value), field);
  }

  constexpr bool operator==(const element_t &other) const {
    return field.eq(value, other.value);
  }

  friend std::ostream &operator<<(std::ostream &os, const FieldElement &val) {
    return os << val.value;
  }
};
