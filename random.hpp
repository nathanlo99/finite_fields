
#pragma once

#include <cstdint>
#include <random>

namespace util {

inline int64_t random_int64_t(const int64_t min, const int64_t max) {
  static std::random_device device;
  static std::mt19937 rng(device());
  std::uniform_int_distribution<std::mt19937::result_type> dist(min, max);
  return dist(rng);
}

}; // namespace util
