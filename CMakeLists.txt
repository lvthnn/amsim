cmake_minimum_required(VERSION 3.16)
project(amsimcpp LANGUAGES CXX VERSION 0.1.0)

# C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_BUILD_TYPE Release)

# Library targets
add_library(amsimcpp STATIC
  src/genome.cc
  src/phenotype.cc
  src/utils.cc
  src/mating.cc
  src/phenobuf.cc
  src/phenoarch.cc
  src/haplobuf.cc
  src/metric.cc
  src/simulation.cc
  src/simulation_builder.cc
)

# Public headers for users of the lib
target_include_directories(amsimcpp
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>)

# Link to BLAS
option(USE_BLAS "Enable code acceleration using BLAS" OFF)

if (APPLE AND USE_BLAS)
  target_compile_definitions(amsimcpp PRIVATE ACCELERATE_NEW_LAPACK USE_BLAS)
  target_link_libraries(amsimcpp PRIVATE "-framework Accelerate")
elseif (UNIX AND USE_BLAS)
  find_package(BLAS REQUIRED)
  target_compile_definitions(amsimcpp PRIVATE USE_BLAS)
  target_link_libraries(PRIVATE BLAS::BLAS)
endif()

# Reasonable warnings (optional)
if (CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
  target_compile_options(amsimcpp PRIVATE -Wall -Wextra -Wpedantic)
elseif (MSVC)
  target_compile_options(amsimcpp PRIVATE /W4)
endif()

# Example program (optional)
add_executable(demo examples/demo_new.cc)
target_link_libraries(demo PRIVATE amsimcpp)

add_executable(test_phenoarch examples/test_phenoarch.cc)
target_link_libraries(test_phenoarch PRIVATE amsimcpp)

add_executable(benchmark_polar_m19937 examples/benchmark_polar_m19937.cc)
target_link_libraries(benchmark_polar_m19937 PRIVATE amsimcpp)

add_executable(unif_int_rng examples/unif_int_rng.cc)
target_link_libraries(unif_int_rng PRIVATE amsimcpp)

# Enable tests later (optional)
enable_testing()

# Install rules (so other projects can find+link amsimcpp)
include(GNUInstallDirs)
install(TARGETS amsimcpp
  EXPORT amsimcppTargets
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(EXPORT amsimcppTargets
  NAMESPACE amsimcpp::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/amsimcpp)
