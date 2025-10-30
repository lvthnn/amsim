#include <cstddef>
#include <cstdint>

namespace amsim::utils {

namespace {
void bitmatrix_swap(
    std::uint64_t matrix[], std::size_t width, std::uint64_t mask) {
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
}  // namespace

void bitmatrix_transpose(std::uint64_t* matrix) {
  std::size_t swap_width = 64;
  std::uint64_t swap_mask = static_cast<std::uint64_t>(-1);
  while (swap_width != 1) {
    swap_width >>= 1;
    swap_mask = swap_mask ^ (swap_mask >> swap_width);
    bitmatrix_swap(matrix, swap_width, swap_mask);
  }
}

// function to generate UNLINKED recombination map with constant MAFs

// function to generate UNLINKED recombination map with random MAFs

// function to generate LINKED recombination map with otherwise constant MAFs

// function to generate LINKED recombination map with random MAFs

// function to generate phenotype effect vector with uniform effect sizes

// function to generate phenotype effect vector with random effect sizes
}  // namespace amsim::utils
