#include <vector>

#if defined(__APPLE__)
  #include <Accelerate/Accelerate.h>
#else
  #include <cblas.h>
  #include <lapacke.h>
#endif

#include <amsim/stats.h>
#include <amsim/utils.h>

int main() {
  std::size_t n_elem = 512000;
  std::vector<double> yo(n_elem, 5.0);
  std::vector<double> ones(n_elem, 1.0);
  std::cout << n_elem << " elements\n";
  amsim::utils::time_step("non-blas sum", [&](){ amsim::stats::sum(n_elem, yo.data(), 1); });
  amsim::utils::time_step("blas sum", [&]() { cblas_ddot(n_elem, yo.data(), 1, ones.data(), 1); });
  amsim::utils::time_step("cor", [&]() { amsim::stats::cor(n_elem, yo.data(), 1, ones.data(), 1); });

  // test matrix correlation
  std::vector<double> X{1, 2, 3, 1, 2, 3, 1, 2, 3};
  std::vector<double> Y{1, 3, 2, 1, 3, 2, 1, 3, 2};
	std::vector<double> R(9);

	amsim::utils::time_step("matrix cor", [&]() {
		amsim::stats::cor(3, 3, 3, X.data(), 3, Y.data(), 3, R.data(), 3);
	});

	for (std::size_t i = 0; i < 3; ++i) {
		for (std::size_t j = 0; j < 3; ++j)
			std::cout << R[3 * i + j] << "\t";
		std::cout << "\n";
	}
}
