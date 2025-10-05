#include <cstddef>
#include <cstdint>
#include <stdexcept>

namespace amsim::utils {
  namespace {
    void bitmatrix_swap(std::uint64_t matrix[], std::size_t width,
                        std::uint64_t mask) {
      std::size_t inner, outer;
      for (outer = 0; outer < 64 / (width * 2); outer++) {
        for (inner = 0; inner < width; inner++) {
          std::uint64_t* x = &matrix[(inner) + (outer * width * 2)];
          std::uint64_t* y = &matrix[(inner + width) + (outer * width * 2)];
          *x = ((*y << width) & mask) ^ *x;
          *y = ((*x & mask) >> width) ^ *y;
          *x = ((*y << width) & mask) ^ *x;
        }
      }
    }
  }

  void bitmatrix_transpose(std::uint64_t* matrix) {
    std::size_t swap_width = 64;
    std::uint64_t swap_mask = ((std::uint64_t)(-1));
    while (swap_width != 1) {
      swap_width >>= 1;
      swap_mask = swap_mask ^ (swap_mask >> swap_width);
      bitmatrix_swap(matrix, swap_width, swap_mask);
    }
  }

  double ddot(std::size_t n, double* x, double* y) {
    double s = 0.0;
    for (std::size_t k = 0; k < n; ++k) s += x[k] * y[k];
    return s;
  }
}
