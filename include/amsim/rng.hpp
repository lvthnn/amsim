#pragma once
#include <array>
#include <cstdint>
#include <chrono>
#include <type_traits>
#include <cmath>
#if __cpp_lib_bitops
  #include <bit>
#endif

namespace amsim::rng {

	inline uint64_t auto_seed(uint64_t seed) {
		if (seed != 0) return seed;
		return static_cast<uint64_t>(
			std::chrono::high_resolution_clock::now().time_since_epoch().count()
		);
	}

	struct Xoshiro256ss {
		std::uint64_t s[4]{};

		static inline uint64_t rotl(uint64_t x, int k) noexcept {
    #if __cpp_lib_bitops
      return std::rotl(x, k);
    #else
      return (x << k) | (x >> (64 - k));
    #endif
		}

		inline std::uint64_t next() noexcept {
			const std::uint64_t result = rotl(s[1] * 5, 7) * 9;
			const std::uint64_t t = s[1] << 17;
			s[2] ^= s[0]; s[3] ^= s[1]; s[1] ^= s[2];
			s[0] ^= s[3]; s[2] ^= t; s[3] = rotl(s[3], 45);
			return result;
		}
	};

	namespace detail {
		inline std::uint64_t splitmix64_step(std::uint64_t &x) noexcept {
			std::uint64_t z = (x += 0x9e3779b97f4a7c15ULL);
			z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
			z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
			return z ^ (z >> 31);
		}
	}

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

	template<int BITS>
	using ThrT = std::conditional_t<BITS==8,  std::uint8_t,
							std::conditional_t<BITS==16, std::uint16_t, std::uint64_t>>;

	inline std::uint64_t lowbits_mask(unsigned k) noexcept {
		if (k == 0)  return 0ULL;
		if (k >= 64) return ~0ULL;
		return (1ULL << k) - 1ULL;
	}

	template<int BITS>
	inline ThrT<BITS> prob_to_thr(double p) noexcept {
		static_assert(BITS==8 || BITS==16 || BITS==64, "BITS must be 8,16,64");
		if (p <= 0.0) return ThrT<BITS>(0);
		if (p >= 1.0) return ThrT<BITS>(~ThrT<BITS>(0));
		long double v = std::ldexp((long double)p, BITS);
		std::uint64_t t = static_cast<std::uint64_t>(v);
		if constexpr (BITS < 64) {
			const std::uint64_t cap = (1ULL << BITS) - 1ULL;
			if (t > cap) t = cap;
		}
		return ThrT<BITS>(t);
	}

	// word generator classes

	template<int BITS = 16>
	struct BernoulliWordConst {
		static_assert(BITS==8 || BITS==16 || BITS==64);
		using T = ThrT<BITS>;

		explicit BernoulliWordConst(std::uint64_t seed = 0x0123456789abcdefULL)
		: Tj(0), rng_(seed_xoshiro(auto_seed(seed))) {}

		inline void set_prob(double p) noexcept {
			if (p < 0.0) p = 0.0; if (p > 1.0) p = 1.0;
			Tj = prob_to_thr<BITS>(p);
		}

		inline void set_probs(const std::array<double,64>& f) noexcept {
			set_prob(f[0]); // constant-p semantics
		}

		inline void reseed(std::uint64_t seed) noexcept { rng_ = seed_xoshiro(seed); }
		inline bool coinflip() noexcept { return rng_.next() & 1ULL; }

		inline std::uint64_t sample(unsigned valid_bits = 64) noexcept {
			std::uint64_t w = 0;
			if constexpr (BITS == 8) {
				for (int j = 0; j < 64; ) {
					std::uint64_t r = rng_.next();
					for (int k = 0; k < 8 && j < 64; ++k, ++j) {
						std::uint8_t rv = static_cast<std::uint8_t>(r >> 56); r <<= 8;
						w |= (std::uint64_t)-(rv < Tj) & (1ULL << j);
					}
				}
			} else if constexpr (BITS == 16) {
				for (int j = 0; j < 64; ) {
					std::uint64_t r = rng_.next();
					for (int k = 0; k < 4 && j < 64; ++k, ++j) {
						std::uint16_t rv = static_cast<std::uint16_t>(r >> 48); r <<= 16;
						w |= (std::uint64_t)-(rv < Tj) & (1ULL << j);
					}
				}
			} else { // 64
				for (int j = 0; j < 64; ++j) {
					std::uint64_t rv = rng_.next();
					w |= (std::uint64_t)-(rv < Tj) & (1ULL << j);
				}
			}
			if (valid_bits < 64) w &= lowbits_mask(valid_bits);
			return w;
		}

	private:
		T Tj;
		Xoshiro256ss rng_;
	};

	template<int BITS = 16>
	struct BernoulliWordVar {
		static_assert(BITS==8 || BITS==16 || BITS==64);
		using T = ThrT<BITS>;

		explicit BernoulliWordVar(std::uint64_t seed = 0x0123456789abcdefULL)
		: rng_(seed_xoshiro(auto_seed(seed))) { Tj.fill(T(0)); }

		inline void set_prob(double p) noexcept {
			if (p < 0.0) p = 0.0; if (p > 1.0) p = 1.0;
			const T t = prob_to_thr<BITS>(p);
			Tj.fill(t);
		}
		inline void set_probs(const std::array<double,64>& f) noexcept {
			for (int j=0; j<64; ++j) {
				double p = f[j]; if (p < 0.0) p = 0.0; if (p > 1.0) p = 1.0;
				Tj[j] = prob_to_thr<BITS>(p);
			}
		}

		inline void reseed(std::uint64_t seed) noexcept { rng_ = seed_xoshiro(seed); }
		inline bool coinflip() noexcept { return rng_.next() & 1ULL; }

		inline std::uint64_t sample(unsigned valid_bits = 64) noexcept {
			std::uint64_t w = 0;
			if constexpr (BITS == 8) {
				for (int j = 0; j < 64; ) {
					std::uint64_t r = rng_.next();
					for (int k = 0; k < 8 && j < 64; ++k, ++j) {
						std::uint8_t rv = static_cast<std::uint8_t>(r >> 56); r <<= 8;
						w |= (std::uint64_t)-(rv < Tj[j]) & (1ULL << j);
					}
				}
			} else if constexpr (BITS == 16) {
				for (int j = 0; j < 64; ) {
					std::uint64_t r = rng_.next();
					for (int k = 0; k < 4 && j < 64; ++k, ++j) {
						std::uint16_t rv = static_cast<std::uint16_t>(r >> 48); r <<= 16;
						w |= (std::uint64_t)-(rv < Tj[j]) & (1ULL << j);
					}
				}
			} else { // 64
				for (int j = 0; j < 64; ++j) {
					std::uint64_t rv = rng_.next();
					w |= (std::uint64_t)-(rv < Tj[j]) & (1ULL << j);
				}
			}
			if (valid_bits < 64) w &= lowbits_mask(valid_bits);
			return w;
		}

	private:
		std::array<T,64> Tj{};
		Xoshiro256ss rng_;
	};

  using BW16 = BernoulliWordVar<16>;
}
