// This file is part of amsim, copyright (C) 2025-2026 Kári Hlynsson.
//
// This program is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your option)
// any later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
// more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#pragma once

#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <type_traits>
#if __cpp_lib_bitops
#include <bit>
#endif

namespace amsim::rng {

inline uint64_t auto_seed(std::optional<uint64_t> seed) {
  if (seed.has_value()) return seed.value();
  return static_cast<uint64_t>(
      std::chrono::high_resolution_clock::now().time_since_epoch().count());
}

class Xoshiro256ss {
 public:
  Xoshiro256ss(const Xoshiro256ss&) = delete;
  Xoshiro256ss& operator=(const Xoshiro256ss&) = delete;

  static Xoshiro256ss& get_instance() {
    thread_local static Xoshiro256ss rng;
    return rng;
  }

  void set_seed(std::uint64_t seed) noexcept {
    std::uint64_t x = seed ? seed : 0x9e3779b97f4a7c15ULL;
    s_[0] = splitMix64Step(x);
    s_[1] = splitMix64Step(x);
    s_[2] = splitMix64Step(x);
    s_[3] = splitMix64Step(x);
    if ((s_[0] | s_[1] | s_[2] | s_[3]) == 0) s_[0] = 1;
  }

  std::uint64_t next() noexcept {
    const std::uint64_t result = rotl(s_[1] * 5, 7) * 9;
    const std::uint64_t t = s_[1] << 17;
    s_[2] ^= s_[0];
    s_[3] ^= s_[1];
    s_[1] ^= s_[2];
    s_[0] ^= s_[3];
    s_[2] ^= t;
    s_[3] = rotl(s_[3], 45);
    return result;
  }

 private:
  explicit Xoshiro256ss() = default;

  std::uint64_t s_[4]{};  ///< Generator state

  static uint64_t rotl(uint64_t x, int k) noexcept {
#if __cpp_lib_bitops
    return std::rotl(x, k);
#else
    return (x << k) | (x >> (64 - k));
#endif
  }

  static std::uint64_t splitMix64Step(std::uint64_t& x) noexcept {
    std::uint64_t z = (x += 0x9e3779b97f4a7c15ULL);
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
  }
};

template <int BITS>
using ThrT = std::conditional_t<
    BITS == 8,
    std::uint8_t,
    std::conditional_t<BITS == 16, std::uint16_t, std::uint64_t>>;

inline std::uint64_t lowbits_mask(unsigned k) noexcept {
  if (k == 0) return 0ULL;
  if (k >= 64) return ~0ULL;
  return (1ULL << k) - 1ULL;
}

template <int BITS>
inline ThrT<BITS> prob_to_thr(double p) noexcept {
  static_assert(
      BITS == 8 || BITS == 16 || BITS == 64, "BITS must be 8, 16, 64");
  if (p <= 0.0) return ThrT<BITS>(0);
  if (p >= 1.0) return ThrT<BITS>(~ThrT<BITS>(0));
  long double v = std::ldexp(static_cast<long double>(p), BITS);
  auto t = static_cast<std::uint64_t>(v);
  if constexpr (BITS < 64) {
    const std::uint64_t cap = (1ULL << BITS) - 1ULL;
    t = std::min(t, cap);
  }
  return ThrT<BITS>(t);
}

template <int BITS = 16>
struct BernoulliWord {
  static_assert(BITS == 8 || BITS == 16 || BITS == 64);
  using T = ThrT<BITS>;  ///< Threshold type

  BernoulliWord() = default;

  void set_prob(double p) noexcept {
    p = std::max(0.0, p);
    p = std::min(p, 1.0);
    T thresh = prob_to_thr<BITS>(p);
    for (std::size_t j = 0; j < 64; ++j) tj_[j] = thresh;
  }

  void set_probs(const double* ptr, std::size_t valid = 64) noexcept {
    double p;
    for (std::size_t j = 0; j < 64; ++j) {
      p = (j < valid) ? *(ptr + j) : 0;
      p = std::max(0.0, p);
      p = std::min(p, 1.0);
      T thresh = prob_to_thr<BITS>(p);
      tj_[j] = thresh;
    }
  }

  bool coinflip() noexcept {
    return (Xoshiro256ss::get_instance().next() & 1ULL) != 0U;
  }

  std::uint64_t sample(unsigned valid_bits = 64) noexcept {
    std::uint64_t w = 0;
    if constexpr (BITS == 8) {
      for (int j = 0; j < 64;) {
        std::uint64_t r = Xoshiro256ss::get_instance().next();
        for (int k = 0; k < 8 && j < 64; ++k, ++j) {
          auto rv = static_cast<std::uint8_t>(r >> 56);
          r <<= 8;
          w |= static_cast<std::uint64_t>(-(rv < tj_[j])) & (1ULL << j);
        }
      }
    } else if constexpr (BITS == 16) {
      for (int j = 0; j < 64;) {
        std::uint64_t r = Xoshiro256ss::get_instance().next();
        for (int k = 0; k < 4 && j < 64; ++k, ++j) {
          auto rv = static_cast<std::uint16_t>(r >> 48);
          r <<= 16;
          w |= static_cast<std::uint64_t>(-(rv < tj_[j])) & (1ULL << j);
        }
      }
    } else {
      for (int j = 0; j < 64; ++j) {
        std::uint64_t rv = Xoshiro256ss::get_instance().next();
        w |= static_cast<std::uint64_t>(-(rv < tj_[j])) & (1ULL << j);
      }
    }
    if (valid_bits < 64) w &= lowbits_mask(valid_bits);
    return w;
  }

 private:
  std::array<T, 64> tj_;  ///< Probability threshold
};

using BW16 = BernoulliWord<16>;

inline double u01_53(const uint64_t x) noexcept {
  return ((x >> 11) + 0.5) * (1.0 / 9007199254740992.0);
}

struct NormalPolar {
  NormalPolar() = default;

  static std::array<double, 2> two() noexcept {
    double u;
    double v;
    double s;
    do {
      u = (2.0 * u01_53(Xoshiro256ss::get_instance().next())) - 1.0;
      v = (2.0 * u01_53(Xoshiro256ss::get_instance().next())) - 1.0;
      s = (u * u) + (v * v);
    } while (s >= 1.0 || s == 0.0);
    const double m = std::sqrt(-2.0 * std::log(s) / s);
    return {u * m, v * m};
  }

  static double single() noexcept { return two()[0]; }

  static void fill(double* out, std::size_t n) noexcept {
    std::size_t i = 0;
    for (; i + 1 < n; i += 2) {
      auto z = two();
      out[i] = z[0];
      out[i + 1] = z[1];
    }
    if (i < n) out[i] = two()[0];
  }

  static std::vector<double> sample(std::size_t n_elem) noexcept {
    std::vector<double> out(n_elem);
    fill(out.data(), n_elem);
    return out;
  }
};

struct UniformRange {
  UniformRange() = default;

  static double sample(double a = 1.0) noexcept {
    return a * u01_53(Xoshiro256ss::get_instance().next());
  }

  static void fill(double* out, std::size_t n, double a = 1.0) noexcept {
    for (std::size_t i = 0; i < n; ++i)
      out[i] = a * u01_53(Xoshiro256ss::get_instance().next());
  }
};

struct UniformIntRange {
  UniformIntRange() = default;

  static std::size_t sample(const std::size_t hi) {
    if (hi == 0) return 0;

    const std::size_t thresh = UINT64_MAX - (UINT64_MAX % hi);
    std::size_t x;

    do {
      x = Xoshiro256ss::get_instance().next();
    } while (x >= thresh);

    return (x % hi);
  }

  static std::size_t sample(const std::size_t lo, const std::size_t hi) {
    if (hi < lo) throw std::runtime_error("[hi, lo) must be non-empty");
    if (hi == lo) return lo;

    const std::size_t thresh = UINT64_MAX - (UINT64_MAX % (hi - lo));
    std::size_t x;

    do {
      x = Xoshiro256ss::get_instance().next();
    } while (x >= thresh);

    return (x % (hi - lo)) + lo;
  }
};

inline void set_seed(std::uint64_t seed) {
  Xoshiro256ss::get_instance().set_seed(seed);
}

}  // namespace amsim::rng
