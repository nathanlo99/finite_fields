
#pragma once

#include "error.hpp"

#include <string>
#include <vector>

template <class Value = int64_t> struct Polynomial {
  using value_t = Value;

  char variable;
  std::vector<value_t> coeffs;

public:
  // TODO: Check for Ring::zero() instead of 0
  Polynomial(const char variable, const std::vector<value_t> &coeffs)
      : variable(variable), coeffs(coeffs) {
    // Enforce having at least the constant coefficient
    if (coeffs.empty())
      this->coeffs = {0};

    // Remove leading zeroes
    while (this->coeffs.size() > 1 && this->coeffs.back() == 0)
      this->coeffs.pop_back();
  }

  inline value_t operator[](const size_t idx) const {
    // TODO: Make Value a Ring and return Ring::zero() instead of 0
    return idx < coeffs.size() ? coeffs[idx] : 0;
  }

  Polynomial operator-() const {
    const size_t degree_plus_one = coeffs.size();
    std::vector<value_t> result_coeffs(degree_plus_one);
    // TOOD: Parallelize this
    for (size_t i = 0; i < degree_plus_one; ++i) {
      result_coeffs[i] = -coeffs[i];
    }
    return Polynomial(variable, result_coeffs);
  }

  friend Polynomial operator+(const Polynomial &a, const Polynomial &b) {
    if (a.variable != b.variable)
      throw math_error()
          << "Cannot add polynomials with two different variables '"
          << a.variable << "' and '" << b.variable << "'";
    const size_t degree_plus_one = std::max(a.coeffs.size(), b.coeffs.size());
    std::vector<value_t> result_coeffs(degree_plus_one);
    // TODO: Parallelize this
    for (size_t i = 0; i < degree_plus_one; ++i) {
      result_coeffs[i] = a[i] + b[i];
    }
    return Polynomial(a.variable, result_coeffs);
  }

  friend Polynomial operator-(const Polynomial &a, const Polynomial &b) {
    return a + (-b);
  }

  friend Polynomial operator*(const value_t k, const Polynomial &p) {
    const size_t degree_plus_one = p.coeffs.size();
    std::vector<value_t> result_coeffs(degree_plus_one);
    // TODO: Parallelize this
    for (size_t i = 0; i < degree_plus_one; ++i) {
      result_coeffs[i] = k * p.coeffs[i];
    }
    return Polynomial(p.variable, result_coeffs);
  }

  friend Polynomial operator*(const Polynomial &p, const value_t k) {
    return k * p;
  }

  friend Polynomial operator*(const Polynomial &a, const Polynomial &b) {
    if (a.variable != b.variable)
      throw math_error()
          << "Cannot multiply polynomials with two different variables '"
          << a.variable << "' and '" << b.variable << "'";
    // TODO: Write a sub-quadratic multiplication algorithm (maybe FFT?)
    const size_t a_degree = a.coeffs.size() - 1;
    const size_t b_degree = b.coeffs.size() - 1;
    const size_t result_degree = a_degree + b_degree;
    // TODO: Ring::zero()
    std::vector<value_t> result_coeffs(result_degree + 1, 0);
    for (size_t i = 0; i <= a_degree; ++i) {
      for (size_t j = 0; j <= b_degree; ++j) {
        result_coeffs[i + j] += a[i] * b[j];
      }
    }
    return Polynomial(a.variable, result_coeffs);
  }

  Polynomial &operator*=(const value_t k) { return *this = *this * k; }

  Polynomial &operator*=(const Polynomial &other) {
    return *this = *this * other;
  }

  Polynomial operator^(uint64_t exp) const {
    Polynomial result(variable, {1}), pow = *this;
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
    // TODO: Make this more sophisticated: skip zero entries
    if (p.coeffs.size() == 1 && p[0] == 0)
      return os << "0";
    bool first = true;
    for (int d = static_cast<int>(p.coeffs.size()) - 1; d >= 0; --d) {
      if (p[d] == 0)
        continue;
      if (first) {
        first = false;
        if (p.coeffs[d] < 0)
          os << "-";
      } else if (p.coeffs[d] < 0) {
        os << " - ";
      } else {
        os << " + ";
      }

      const value_t coeff = std::abs(p.coeffs[d]);
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
};
