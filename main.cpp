
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>

#include "elliptic_curve.hpp"
#include "number_theory.hpp"
#include "prime_field.hpp"
#include "timing.hpp"

int main() {
  PrimeField F = PrimeField(1000003LL);
  const auto E = EllipticCurve(F, 2, 3);
  std::cout << E << std::endl;

  const auto start_ms = get_time_ms();
  const auto points = E.points();
  std::cout << points.size() << " points" << std::endl;
  const auto end_ms = get_time_ms();
  std::cout << "Done in " << (end_ms - start_ms) / 1000.0 << " seconds"
            << std::endl;

  const auto point = points[points.size() / 2];
  std::cout << point << std::endl;
  std::cout << point + point << std::endl;
  std::cout << point + point + point << std::endl;

  // for (long long p = 3; p < 1'000'000; ++p) {
  //   if (!NumberTheory::is_prime_slow(p))
  //     continue;
  //   const long long expected_num_quadratic_residues = (p + 1) / 2;
  //   long long num_quadratic_residues = 0;
  //   for (long long i = 0; i < p; ++i) {
  //     if (!nt::is_quadratic_residue(i, p))
  //       continue;
  //     num_quadratic_residues++;
  //     const long long sqrt = nt::sqrt_mod(i, p);
  //     const long long prod = nt::mul_mod(sqrt, sqrt, p);
  //     if (prod != i) {
  //       std::cout << "sqrt(" << i << ") failed" << std::endl;
  //       exit(1);
  //     }
  //   }
  //   if (num_quadratic_residues != expected_num_quadratic_residues) {
  //     throw std::runtime_error(
  //         "Wrong number of quadratic residues for p = " + std::to_string(p) +
  //         ", expected " + std::to_string(expected_num_quadratic_residues) +
  //         ", got " + std::to_string(num_quadratic_residues));
  //   }
  //   std::cout << "Done testing " << p << std::endl;
  // }
}
