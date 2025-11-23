#ifndef AMSIMCPP_RNG_H
#define AMSIMCPP_RNG_H

#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <type_traits>

#if __cpp_lib_bitops
#include <bit>
#endif

namespace amsim::rng {

inline uint64_t auto_seed(uint64_t seed) {
  if (seed != 0) return seed;
  return static_cast<uint64_t>(
      std::chrono::high_resolution_clock::now().time_since_epoch().count());
}

struct Xoshiro256ss {
  std::uint64_t s[4]{};

  static uint64_t rotl(uint64_t x, int k) noexcept {
#if __cpp_lib_bitops
    return std::rotl(x, k);
#else
    return (x << k) | (x >> (64 - k));
#endif
  }

  std::uint64_t next() noexcept {
    const std::uint64_t result = rotl(s[1] * 5, 7) * 9;
    const std::uint64_t t = s[1] << 17;
    s[2] ^= s[0];
    s[3] ^= s[1];
    s[1] ^= s[2];
    s[0] ^= s[3];
    s[2] ^= t;
    s[3] = rotl(s[3], 45);
    return result;
  }
};

namespace detail {
inline std::uint64_t splitmix64_step(std::uint64_t& x) noexcept {
  std::uint64_t z = (x += 0x9e3779b97f4a7c15ULL);
  z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
  z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
  return z ^ (z >> 31);
}
}  // namespace detail

inline Xoshiro256ss seed_xoshiro(std::uint64_t seed) noexcept {
  Xoshiro256ss g{};
  std::uint64_t x = seed ? seed : 0x9e3779b97f4a7c15ULL;
  g.s[0] = detail::splitmix64_step(x);
  g.s[1] = detail::splitmix64_step(x);
  g.s[2] = detail::splitmix64_step(x);
  g.s[3] = detail::splitmix64_step(x);
  if ((g.s[0] | g.s[1] | g.s[2] | g.s[3]) == 0) g.s[0] = 1;
  return g;
}

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

// word generator classes

template <int BITS = 16>
struct BernoulliWordConst {
  static_assert(BITS == 8 || BITS == 16 || BITS == 64);
  using T = ThrT<BITS>;

  explicit BernoulliWordConst(const Xoshiro256ss& rng) : tj_(0), rng_(rng) {}

  void set_prob(double p) noexcept {
    p = std::max(0.0, p);
    p = std::min(p, 1.0);
    tj_ = prob_to_thr<BITS>(p);
  }

  void set_probs(const std::array<double, 64>& f) noexcept {
    set_prob(f[0]);  // constant-p semantics
  }

  void reseed(std::uint64_t seed) noexcept { rng_ = seed_xoshiro(seed); }
  bool coinflip() noexcept { return rng_.next() & 1ULL; }

  std::uint64_t sample(unsigned valid_bits = 64) noexcept {
    std::uint64_t w = 0;
    if constexpr (BITS == 8) {
      for (int j = 0; j < 64;) {
        std::uint64_t r = rng_.next();
        for (int k = 0; k < 8 && j < 64; ++k, ++j) {
          auto rv = static_cast<std::uint8_t>(r >> 56);
          r <<= 8;
          w |= static_cast<std::uint64_t>(-(rv < tj_)) & (1ULL << j);
        }
      }
    } else if constexpr (BITS == 16) {
      for (int j = 0; j < 64;) {
        std::uint64_t r = rng_.next();
        for (int k = 0; k < 4 && j < 64; ++k, ++j) {
          auto rv = static_cast<std::uint16_t>(r >> 48);
          r <<= 16;
          w |= static_cast<std::uint64_t>(-(rv < tj_)) & (1ULL << j);
        }
      }
    } else {  // 64
      for (int j = 0; j < 64; ++j) {
        std::uint64_t rv = rng_.next();
        w |= static_cast<std::uint64_t>(-(rv < tj_)) & (1ULL << j);
      }
    }
    if (valid_bits < 64) w &= lowbits_mask(valid_bits);
    return w;
  }

 private:
  T tj_;
  Xoshiro256ss rng_;
};

template <int BITS = 16>
struct BernoulliWordVar {
  static_assert(BITS == 8 || BITS == 16 || BITS == 64);
  using T = ThrT<BITS>;

  explicit BernoulliWordVar(const Xoshiro256ss& rng) : rng_(rng) {
    tj_.fill(T(0));
  }

  void set_prob(double p) noexcept {
    p = std::max(0.0, p);
    p = std::min(p, 1.0);
    const T t = prob_to_thr<BITS>(p);
    tj_.fill(t);
  }

  void set_probs(const std::array<double, 64>& f) noexcept {
    for (int j = 0; j < 64; ++j) {
      double p = f[j];
      p = std::max(0.0, p);
      p = std::min(p, 1.0);
      tj_[j] = prob_to_thr<BITS>(p);
    }
  }

  void reseed(std::uint64_t seed) noexcept { rng_ = seed_xoshiro(seed); }
  bool coinflip() noexcept { return rng_.next() & 1ULL; }

  std::uint64_t sample(unsigned valid_bits = 64) noexcept {
    std::uint64_t w = 0;
    if constexpr (BITS == 8) {
      for (int j = 0; j < 64;) {
        std::uint64_t r = rng_.next();
        for (int k = 0; k < 8 && j < 64; ++k, ++j) {
          auto rv = static_cast<std::uint8_t>(r >> 56);
          r <<= 8;
          w |= static_cast<std::uint64_t>(-(rv < tj_[j])) & (1ULL << j);
        }
      }
    } else if constexpr (BITS == 16) {
      for (int j = 0; j < 64;) {
        std::uint64_t r = rng_.next();
        for (int k = 0; k < 4 && j < 64; ++k, ++j) {
          auto rv = static_cast<std::uint16_t>(r >> 48);
          r <<= 16;
          w |= static_cast<std::uint64_t>(-(rv < tj_[j])) & (1ULL << j);
        }
      }
    } else {  // 64
      for (int j = 0; j < 64; ++j) {
        std::uint64_t rv = rng_.next();
        w |= static_cast<std::uint64_t>(-(rv < tj_[j])) & (1ULL << j);
      }
    }
    if (valid_bits < 64) w &= lowbits_mask(valid_bits);
    return w;
  }

 private:
  std::array<T, 64> tj_{};
  Xoshiro256ss rng_;
};

using BW16 = BernoulliWordConst<16>;

// polar method for generating iid standard normal variables

inline double u01_53(const uint64_t x) noexcept {
  return ((x >> 11) + 0.5) * (1.0 / 9007199254740992.0);
}

struct NormalPolar {
  Xoshiro256ss rng;

  explicit NormalPolar(const Xoshiro256ss& rng_) : rng(rng_) {}

  void reseed(uint64_t seed) noexcept { rng = seed_xoshiro(seed); }

  std::array<double, 2> two() noexcept {
    double u;
    double v;
    double s;
    do {
      u = 2.0 * u01_53(rng.next()) - 1.0;
      v = 2.0 * u01_53(rng.next()) - 1.0;
      s = u * u + v * v;
    } while (s >= 1.0 || s == 0.0);
    const double m = std::sqrt(-2.0 * std::log(s) / s);
    return {u * m, v * m};
  }

  void fill(double* out, std::size_t n) noexcept {
    std::size_t i = 0;
    for (; i + 1 < n; i += 2) {
      auto z = two();
      out[i] = z[0];
      out[i + 1] = z[1];
    }
    if (i < n) out[i] = two()[0];
  }

  template <std::size_t N>
  std::array<double, N> batch() noexcept {
    static_assert(N % 2 == 0, "N must be even");
    std::array<double, N> a{};
    for (std::size_t i = 0; i < N; i += 2) {
      auto z = two();
      a[i] = z[0];
      a[i + 1] = z[1];
    }
    return a;
  }
};

// @TODO: Convert structs `UniformRange` and `UniformIntRange` into
//        generalised intervals [a, b].
struct UniformRange {
  Xoshiro256ss rng;

  explicit UniformRange(const Xoshiro256ss& rng_) : rng(rng_) {}

  void reseed(uint64_t seed) noexcept { rng = seed_xoshiro(seed); }

  double sample(double a) noexcept { return a * u01_53(rng.next()); }

  void fill(double* out, std::size_t n, double a) noexcept {
    for (std::size_t i = 0; i < n; ++i) out[i] = a * u01_53(rng.next());
  }
};

struct UniformIntRange {
  Xoshiro256ss rng;

  explicit UniformIntRange(const Xoshiro256ss& rng_) : rng(rng_) {}

  std::size_t sample(const std::size_t n) {
    if (n == 0) return 0;

    const std::size_t thresh = UINT64_MAX - (UINT64_MAX % n);
    std::size_t x;

    do {
      x = rng.next();
    } while (x >= thresh);

    return x % n;
  }
};

}  // namespace amsim::rng

#endif  // AMSIMCPP_RNG_H
