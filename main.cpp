
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

  timeit("generate_elliptic_curve_points", [&]() {
    const auto points = E.points();
    std::cout << points.size() << " points" << std::endl;

    const auto point = points[points.size() / 2];
    std::cout << point << std::endl;
    std::cout << point + point << std::endl;
    std::cout << point + point + point << std::endl;
  });

  timeit("check_square_roots", []() {
    for (long long p = 3; p < 5'000; ++p) {
      if (!NumberTheory::is_prime(p))
        continue;
      const long long expected_num_quadratic_residues = (p + 1) / 2;
      long long num_quadratic_residues = 0;
      for (long long i = 0; i < p; ++i) {
        if (!nt::is_quadratic_residue(i, p))
          continue;
        num_quadratic_residues++;
        const long long sqrt = nt::sqrt_mod(i, p);
        const long long prod = nt::mul_mod(sqrt, sqrt, p);
        if (prod != i)
          throw std::runtime_error("sqrt(" + std::to_string(i) + " failed");
      }
      if (num_quadratic_residues != expected_num_quadratic_residues) {
        throw std::runtime_error(
            "Wrong number of quadratic residues for p = " + std::to_string(p) +
            ", expected " + std::to_string(expected_num_quadratic_residues) +
            ", got " + std::to_string(num_quadratic_residues));
      }
    }
  });

  timeit("prime_count", []() {
    long long num_primes = 0, max_num = 1'000'000;
    for (long long p = 0; p < max_num; ++p) {
      if (NumberTheory::is_prime_miller_rabin(p))
        num_primes++;
    }
    std::cout << "There are " << num_primes << " primes under " << max_num
              << std::endl;
  });
}
