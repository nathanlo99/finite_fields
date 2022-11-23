
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>

#include "elliptic_curve.hpp"
#include "fraction.hpp"
#include "number_theory.hpp"
#include "polynomial.hpp"
#include "prime_field.hpp"
#include "rational_field.hpp"
#include "timing.hpp"

int main() {
  timeit("rational_polynomial", [&]() {
    const auto QQ = RationalField<int64_t>();
    const auto x = Polynomial(QQ, 'x');
    std::cout << ((x + 1) ^ 5) << std::endl;
    const auto y = Polynomial(QQ, 'y');
    std::cout << ((y + 2) ^ 8) << std::endl;
    std::cout << (x ^ 2) + 2 * x + 1 << std::endl;

    const auto f = Polynomial(QQ, 'x', {1, 2, 3});
    const auto g = Polynomial(QQ, 'x', {1, 1, 4});

    std::cout << "f = " << f << std::endl;
    std::cout << "g = " << g << std::endl;
    std::cout << "f / g = " << f / g << std::endl;
    std::cout << "f % g = " << f % g << std::endl;

    const auto F5 = PrimeField<int64_t>(5);
    const auto h1 = Polynomial(F5, 'x', {1, 2, 3});
    const auto h2 = Polynomial(F5, 'x', {1, 1, 4});
    std::cout << "h1 = " << h1 << std::endl;
    std::cout << "h2 = " << h2 << std::endl;
    std::cout << "h1 / h2 = " << h1 / h2 << std::endl;
    std::cout << "h1 % h2 = " << h1 % h2 << std::endl;
  });

  timeit("generate_elliptic_curve_points", [&]() {
    PrimeField F = PrimeField(1000003LL);
    const auto E = EllipticCurve(F, 2, 3);
    std::cout << E << std::endl;
    const auto points = E.points();
    std::cout << points.size() << " points" << std::endl;

    const auto point = points[points.size() / 2];
    auto sum = E.infinity();
    for (int i = 0; i < 1'000'000; ++i) {
      const auto product = i * point;
      assert(product == sum);
      sum = sum + point;
    }
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
    long long num_primes = 0, max_num = 10'000'000;
    for (long long p = 2; p < max_num; p = nt::next_prime(p)) {
      num_primes++;
    }
    std::cout << "There are " << num_primes << " primes under " << max_num
              << std::endl;
  });
}
