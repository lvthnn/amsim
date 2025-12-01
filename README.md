# amsimcpp: A high-performance C++ library for forward-time assortative mating simulations

**amsimcpp** is a high-performance C++ library for forward-time population
genetic simulations under assortative mating. It provides a computation core for
generating data across multiple generations with configurable genetic
architectures and mating patterns.

### Features

- Multi-trait modelling with customisable phenogenetic architecture
- Assortative and random mating regimes
- Parallel simulation execution with C++ native multithreading
- Efficient linear algebra operations via BLAS/LAPACK

### Language bindings

To run simulations, end users should consult the language-specific bindings.

- **Python:** [amsimpy](./python/) - Install locally via `pip` 
- **R:** [amsimr](./r/) - Install via `devtools::install_github` or locally

*Note: PyPI and CRAN release are planned for future versions.*

These bindings provide interfaces to the C++ simulation core. See the
respective repositories for installation instructions and usage examples.

### Building from source

**amsimcpp** uses CMake and requires BLAS/LAPACK (or the macOS Accelerate
framework, installed by default on newer machines).

```bash
git clone https://github.com/lvthnn/amsimcpp
cd amsimcpp
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
```

We highly recommend compiling from source with `DCMAKE_BUILD_TYPE=Release`, as
simulation performance is significantly improved by the optimisation performed.
However, the `Debug` flag is useful for development to catch segmentation faults
and the like since it builds with sanitisiers enabled. To build the Python
binding `.so` file, specify `-DBUILD_PYTHON=ON` (although this is handled by the
language itself when installing the package).

To run unit tests:
```bash
cd build
ctest --rerun-failed --output-on-failure
```

### Requirements

- C++20 compatible compiler
- CMake ≥ 3.20
- BLAS/LAPACK vendor 

### Examples

See the [examples directory](./examples/) for example programs:

- `demo_simulation.cc` - Basic simulation configuration and execution
- `demo_results.cc` - Working with simulation results
- `demo_logger.cc` - Logging and debugging utilities

Examples are built by default. However, you can also do

```bash
cmake --build build --target demo_simulation
./build/demo_simulation
```

### Installation

Standard CMake package installation is supported:
```bash
cmake --install build --prefix /usr/local
```
Other CMake projects can then find and link against **amsimcpp**:
```
find_package(amsimcpp REQUIRED)
target_link_libraries(your_target PRIVATE amsimcpp::amsimcpp)
```

### License

MIT license - see [LICENSE](./LICENSE) file for details.

### Citation

If you use amsimcpp in your research, please cite:

```bibtex
@software{hlynsson2024_amsimcpp,
  author = {Hlynsson, Kári},
  title = {amsimcpp: A C++ library for forward-time assortative mating simulations},
  year = {2025},
  version = {0.1.0},
  url = {https://github.com/lvthnn/amsimcpp}
}
```

