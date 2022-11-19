
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>

#include "elliptic_curve.hpp"
#include "number_theory.hpp"
#include "prime_field.hpp"

namespace nt = NumberTheory;

int main() {
  PrimeField<> F = PrimeField<>(1000003);
  // std::cout << F << std::endl;

  // std::cout << "       ";
  // for (int j = 0; j < F.p; ++j) {
  //   std::cout << std::setw(3) << j << " ";
  // }
  // std::cout << std::endl;
  // std::cout << "     +-";
  // for (int j = 0; j < F.p; ++j) {
  //   std::cout << "----";
  // }
  // std::cout << std::endl;

  // for (int i = 0; i < F.p; ++i) {
  //   std::cout << std::setw(4) << i << " | ";
  //   for (int j = 0; j < F.p; ++j) {
  //     const auto a = F.value(i);
  //     const auto b = F.value(j);
  //     std::cout << std::setw(3) << a * b << " ";
  //   }
  //   std::cout << std::endl;
  // }

  const auto E = EllipticCurve(F, 2, 3);
  std::cout << E << std::endl;

  int num_points = 0;
  for (long long x = 0; x < F.p; ++x) {
    const auto x_element = F(x);
    const auto fx = E.f(x_element);
    if (!nt::is_quadratic_residue(fx.value, F.p))
      continue;

    const auto y_element = F(nt::sqrt_mod(fx.value, F.p));
    const auto p1 = E.affine_point(x_element, y_element);
    const auto p2 = E.affine_point(x_element, -y_element);
    std::cout << p1 << std::endl;
    num_points++;
    if (p1 != p2) {
      std::cout << p2 << std::endl;
      num_points++;
    }
  }
  std::cout << num_points + 1 << " points" << std::endl;
  // 999705 points

  return 0;

  for (long long p = 2; p < 1000000; ++p) {
    if (!NumberTheory::is_prime_slow(p))
      continue;
    for (long long i = 0; i < p; ++i) {
      if (!nt::is_quadratic_residue(i, p))
        continue;
      const long long sqrt = nt::sqrt_mod(i, p);
      const long long prod = nt::mul_mod(sqrt, sqrt, p);
      if (prod != i) {
        std::cout << "sqrt(" << i << ") failed" << std::endl;
        exit(1);
      }
    }
    std::cout << "Done testing " << p << std::endl;
  }
}
