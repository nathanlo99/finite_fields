
#pragma once

#include "error.hpp"
#include "rational_field.hpp"

#include <string>
#include <vector>

template <class Field = RationalField<int64_t>> struct Polynomial {
  using value_t = typename Field::value_t;
  using element_t = typename Field::element_t;

  const Field &field;
  const element_t zero;
  char variable;
  std::vector<element_t> coeffs;

public:
  Polynomial(const Field &field, const char variable,
             const std::vector<element_t> &coeffs = {})
      : field(field), zero(element_t(field.zero(), field)), variable(variable),
        coeffs(coeffs) {
    if (coeffs.empty()) {
      this->coeffs = {element_t(field.zero(), field),
                      element_t(field.one(), field)};
      return;
    }

    // Remove leading zeroes
    while (this->coeffs.size() > 1 && this->coeffs.back() == zero)
      this->coeffs.pop_back();
  }

  Polynomial(const Field &field, const char variable,
             const std::vector<value_t> &coeffs)
      : field(field), zero(element_t(field.zero(), field)), variable(variable),
        coeffs(coeffs.size(), zero) {
    if (coeffs.empty()) {
      this->coeffs = {element_t(field.zero(), field),
                      element_t(field.one(), field)};
      return;
    }

    for (size_t i = 0; i < coeffs.size(); ++i)
      this->coeffs[i] = element_t(coeffs[i], field);

    // Remove leading zeroes
    while (this->coeffs.size() > 1 && this->coeffs.back() == zero)
      this->coeffs.pop_back();
  }

  Polynomial(const Polynomial &other)
      : field(other.field), zero(other.zero), variable(other.variable),
        coeffs(other.coeffs) {}

  Polynomial &operator=(const Polynomial &other) {
    if (field != other.field)
      throw math_error(
          "Cannot assign Polynomial to Polynomial with different field");
    if (variable != other.variable)
      throw math_error(
          "Cannot assign Polynomial to Polynomial with different variable");
    coeffs = other.coeffs;
    return *this;
  }

  inline element_t operator[](const size_t idx) const {
    return idx < coeffs.size() ? coeffs[idx] : zero;
  }

  Polynomial operator-() const {
    const size_t degree_plus_one = coeffs.size();
    std::vector<element_t> result_coeffs(degree_plus_one, zero);
    // TOOD: Parallelize this
    for (size_t i = 0; i < degree_plus_one; ++i) {
      result_coeffs[i] = -coeffs[i];
    }
    return Polynomial(field, variable, result_coeffs);
  }

  friend Polynomial operator+(const Polynomial &a, const Polynomial &b) {
    if (a.field != b.field)
      throw math_error()
          << "Cannot add polynomials with elements from different fields";
    if (a.variable != b.variable)
      throw math_error()
          << "Cannot add polynomials with two different variables '"
          << a.variable << "' and '" << b.variable << "'";
    const size_t degree_plus_one = std::max(a.coeffs.size(), b.coeffs.size());
    std::vector<element_t> result_coeffs(degree_plus_one, a.zero);
    // TODO: Parallelize this
    for (size_t i = 0; i < degree_plus_one; ++i) {
      result_coeffs[i] = a[i] + b[i];
    }
    return Polynomial(a.field, a.variable, result_coeffs);
  }

  friend Polynomial operator-(const Polynomial &a, const Polynomial &b) {
    return a + (-b);
  }

  friend Polynomial operator*(const value_t k, const Polynomial &p) {
    const size_t degree_plus_one = p.coeffs.size();
    std::vector<element_t> result_coeffs(degree_plus_one, p.zero);
    // TODO: Parallelize this
    for (size_t i = 0; i < degree_plus_one; ++i) {
      result_coeffs[i] = k * p.coeffs[i];
    }
    return Polynomial(p.field, p.variable, result_coeffs);
  }

  friend Polynomial operator*(const Polynomial &p, const value_t k) {
    return k * p;
  }

  friend Polynomial operator*(const Polynomial &a, const Polynomial &b) {
    if (a.field != b.field)
      throw math_error()
          << "Cannot multiply polynomials with elements from different fields";
    if (a.variable != b.variable)
      throw math_error()
          << "Cannot multiply polynomials with two different variables '"
          << a.variable << "' and '" << b.variable << "'";
    // TODO: Write a sub-quadratic multiplication algorithm (maybe FFT?)
    const size_t a_degree = a.coeffs.size() - 1;
    const size_t b_degree = b.coeffs.size() - 1;
    const size_t result_degree = a_degree + b_degree;
    std::vector<element_t> result_coeffs(result_degree + 1, a.zero);
    for (size_t i = 0; i <= a_degree; ++i) {
      for (size_t j = 0; j <= b_degree; ++j) {
        result_coeffs[i + j] += a[i] * b[j];
      }
    }
    return Polynomial(a.field, a.variable, result_coeffs);
  }

  Polynomial &operator*=(const value_t k) { return *this = *this * k; }

  Polynomial &operator*=(const Polynomial &other) {
    return *this = *this * other;
  }

  Polynomial operator^(uint64_t exp) const {
    Polynomial result(field, variable, {1}), pow = *this;
    while (exp > 0) {
      if (exp % 2 == 1)
        result *= pow;
      pow *= pow;
      exp /= 2;
    }
    return result;
  }

  friend constexpr bool operator==(const Polynomial &a, const Polynomial &b) {
    return a.variable == b.variable && a.coeffs == b.coeffs;
  }

  friend std::ostream &operator<<(std::ostream &os, const Polynomial &p) {
    if (p.coeffs.size() == 1 && p.coeffs[0] == p.zero)
      return os << "0";
    bool first = true;
    for (int d = static_cast<int>(p.coeffs.size()) - 1; d >= 0; --d) {
      if (p.coeffs[d] == p.zero)
        continue;
      if (first) {
        first = false;
      } else {
        os << " + ";
      }
      os << p.coeffs[d];
      if (d > 0)
        os << p.variable;
      if (d > 1)
        os << "^" << d;
    }
    return os;
  }
  /*
  friend std::ostream &operator<<(std::ostream &os, const Polynomial &p) {
    if (p.coeffs.size() == 1 && p[0] == 0)
      return os << "0";
    bool first = true;
    for (int d = static_cast<int>(p.coeffs.size()) - 1; d >= 0; --d) {
      if (p[d] == 0)
        continue;
      if (first) {
        first = false;
        if (p.coeffs[d].value < 0)
          os << "-";
      } else if (p.coeffs[d] < 0) {
        os << " - ";
      } else {
        os << " + ";
      }

      const value_t coeff = std::abs(p.coeffs[d].value);
      if (d != 0) {
        if (coeff != 1)
          os << coeff;
        os << p.variable;
      } else {
        os << coeff;
      }
      if (d > 1)
        os << "^" << d;
    }
    return os;
  }
  */
};
