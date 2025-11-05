#ifndef AMSIMCPP_UTILS_H
#define AMSIMCPP_UTILS_H

#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <numeric>

namespace amsim::utils {
void bitmatrix_transpose(std::uint64_t* matrix);

template <typename T>
std::vector<std::size_t> order(const std::vector<T>& v) {
  std::vector<std::size_t> idx(v.size());
  std::iota(idx.begin(), idx.end(), 0);
  std::stable_sort(
      idx.begin(), idx.end(), [&v](std::size_t i0, std::size_t i1) {
        return v[i0] < v[i1];
      });
  return idx;
}

void random_effects(std::size_t n_loci);
void random_mafs(std::size_t n_loci);
void uniform_effects(std::size_t n_loci);
void uniform_mafs(std::size_t n_loci);

template <class F>
void time_step(const char* label, F&& fn) {
  const auto t0 = std::chrono::steady_clock::now();
  fn();
  const auto t1 = std::chrono::steady_clock::now();
  const auto ms =
      std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();
  std::cout << label << " (" << ms << " ms)\n";
}
}  // namespace amsim::utils

#endif  // AMSIMCPP_UTILS_H
