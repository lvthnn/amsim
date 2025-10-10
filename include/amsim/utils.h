#ifndef AMSIMCPP_UTILS_H
#define AMSIMCPP_UTILS_H

#pragma once
#include <cstddef>
#include <cstdint>
#include <iostream>

namespace amsim::utils {
  void bitmatrix_transpose(std::uint64_t* matrix);

  void random_effects(std::size_t n_loci);
  void random_mafs(std::size_t n_loci);
  void uniform_effects(std::size_t n_loci);
  void uniform_mafs(std::size_t n_loci);

  template<class F>
  void time_step(const char *label, F &&fn) {
    const auto t0 = std::chrono::steady_clock::now();
    fn();
    const auto t1 = std::chrono::steady_clock::now();
    const auto ms =
        std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();
    std::cout << label << " (" << ms << " ms)\n";
  }
}

#endif // AMSIMCPP_UTILS_H