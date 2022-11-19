
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

  // const auto E = EllipticCurve(&F, 2, 3);
  // std::cout << E << std::endl;

  // int num_points = 0;
  // for (int x = 0; x < F.p; ++x) {
  //   for (int y = 0; y < F.p; ++y) {
  //     const auto point = E.affine_point(x, y);
  //     if (E.is_point_on_curve(point)) {
  //       std::cout << point << std::endl;
  //       num_points++;
  //     }
  //   }
  // }
  // std::cout << num_points + 1 << " points" << std::endl;

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
