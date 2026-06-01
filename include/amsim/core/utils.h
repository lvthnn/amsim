#pragma once

#include <amsim/core/log.h>

#include <Eigen/Dense>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <numeric>
#include <sys/wait.h>
#include <semaphore>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

namespace amsim::utils {

inline void bitmatrix_swap(
    std::uint64_t matrix[], std::size_t width, std::uint64_t mask) {
  std::size_t inner;
  std::size_t outer;
  for (outer = 0; outer < 64 / (width * 2); ++outer) {
    for (inner = 0; inner < width; ++inner) {
      std::uint64_t* x = &matrix[(inner) + (outer * width * 2)];
      std::uint64_t* y = &matrix[(inner + width) + (outer * width * 2)];
      *x = ((*y << width) & mask) ^ *x;
      *y = ((*x & mask) >> width) ^ *y;
      *x = ((*y << width) & mask) ^ *x;
    }
  }
}

inline void bitmatrix_transpose(std::uint64_t* matrix) {
  std::size_t swap_width = 64;
  auto swap_mask = static_cast<std::uint64_t>(-1);
  while (swap_width != 1) {
    swap_width >>= 1;
    swap_mask = swap_mask ^ (swap_mask >> swap_width);
    bitmatrix_swap(matrix, swap_width, swap_mask);
  }
}

inline std::vector<std::string> vector_suffix(
    const std::vector<std::string>& labels,
    const std::optional<std::string>& suffix = std::nullopt) {
  std::vector<std::string> result(labels.size());
  for (std::size_t i = 0; i < labels.size(); ++i)
    result[i] = labels[i] + suffix.value_or("");
  return result;
}

inline std::vector<std::size_t> order(const Eigen::VectorXd& v) {
  std::vector<std::size_t> idx(v.size());
  std::iota(idx.begin(), idx.end(), 0);
  std::ranges::stable_sort(
      idx, [&v](std::size_t i0, std::size_t i1) { return v(i0) < v(i1); });
  return idx;
}

inline Eigen::MatrixXd standardise(
    const Eigen::MatrixXd& mat, bool population = true) {
  Eigen::RowVectorXd mean = mat.colwise().mean();
  Eigen::RowVectorXd std =
      ((mat.rowwise() - mean).array().square().colwise().sum() /
       (population ? mat.rows() : mat.rows() - 1))
          .sqrt()
          .max(1e-10);
  return ((mat.rowwise() - mean).array().rowwise() / std.array());
}

inline std::counting_semaphore<64> process_semaphore{
    static_cast<std::ptrdiff_t>(
        std::thread::hardware_concurrency() > 0
            ? std::thread::hardware_concurrency()
            : 4)};

inline void system_throttled(const std::string& cmd) {
  process_semaphore.acquire();
  FILE* pipe = popen((cmd + " 2>&1").c_str(), "r");
  std::string output;
  if (pipe) {
    char buffer[512];
    while (fgets(buffer, sizeof(buffer), pipe))
      output += buffer;
  }
  int rc = pclose(pipe);
  process_semaphore.release();
  if (WIFEXITED(rc) && WEXITSTATUS(rc) != 0) {
    Log::error(output);
    throw std::runtime_error(output);
  }
}

inline void check_plink2() {
  int rc = std::system("command -v plink2 >/dev/null 2>&1");
  if (!WIFEXITED(rc) || WEXITSTATUS(rc) != 0)
    throw std::runtime_error("Executable 'plink2' not found in PATH");
}

inline void check_gcta64() {
  int rc = std::system("command -v gcta64 >/dev/null 2>&1");
  if (!WIFEXITED(rc) || WEXITSTATUS(rc) != 0)
    throw std::runtime_error("Executable 'gcta64' not found in PATH");
}

}  // namespace amsim::utils
