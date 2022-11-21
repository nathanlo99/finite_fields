
#pragma once

#include <chrono>
#include <iostream>

inline long long get_time_ns() {
  return std::chrono::duration_cast<std::chrono::nanoseconds>(
             std::chrono::system_clock::now().time_since_epoch())
      .count();
}

inline long long get_time_ms() {
  return std::chrono::duration_cast<std::chrono::milliseconds>(
             std::chrono::system_clock::now().time_since_epoch())
      .count();
}

template <typename Func>
inline void timeit(const std::string_view &msg, const Func &f) {
  const auto start_ms = get_time_ms();
  f();
  const auto end_ms = get_time_ms();
  std::cout << "[" << msg << "] | Done in " << (end_ms - start_ms) / 1000.0
            << " seconds" << std::endl;
}
