# `amsimcpp`: Forward time assortative mating simulations in C++ 
The `amsimcpp` library is written to perform forward-time genomic populations
under assortative mating, with extensible mating regimes and optimisation
methods.

To run `example/demo.cpp`, clone `amsimcpp` and build the package. For
accelerated phenotype scoring using BLAS, run `cmake` with the `-DUSE_BLAS`
option `ON`.
```
git clone https://github.com/lvthnn/amsimcpp.git
cd amsimcpp
mkdir build && cd build
# if you have BLAS:
cmake -S .. -B . -DUSE_BLAS=ON 
# else:
cmake -S .. -B .
cmake --build . --target demo
./demo
```
