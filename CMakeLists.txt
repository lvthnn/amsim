cmake_minimum_required(VERSION 3.20)
project(amsimcpp LANGUAGES CXX VERSION 0.1.0)

# C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Add debugging options depending on CMAKE_BUILD_TYPE
if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    add_compile_options(-fsanitize=address,undefined -fno-omit-frame-pointer)
    add_link_options(-fsanitize=address,undefined -fno-omit-frame-pointer)
endif()

# Register amsim as a header-only library
add_library(amsimcpp INTERFACE)

# Public headers for users of the lib
target_include_directories(amsimcpp
  INTERFACE
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>)

# Use Eigen with BLAS/LAPACK

if (APPLE)
  find_library(ACCELERATE_FRAMEWORK Accelerate REQUIRED)
  target_link_libraries(amsimcpp INTERFACE ${ACCELERATE_FRAMEWORK})
  target_compile_definitions(amsimcpp INTERFACE ACCELERATE_NEW_LAPACK)
  target_compile_definitions(amsimcpp INTERFACE EIGEN_USE_BLAS EIGEN_USE_LAPACK)
else()
  include_directories("/usr/include/openblas")
  find_package(BLAS REQUIRED)
  find_package(LAPACK REQUIRED)
  target_link_libraries(amsimcpp INTERFACE BLAS::BLAS LAPACK::LAPACK)
  target_compile_definitions(amsimcpp INTERFACE EIGEN_USE_BLAS EIGEN_USE_LAPACK)
endif()

# Use highfive for HDF5 support

find_package(HDF5 REQUIRED COMPONENTS C)
target_link_libraries(amsimcpp INTERFACE HDF5::HDF5)
target_compile_definitions(amsimcpp INTERFACE HIGHFIVE_USE_EIGEN)

# Python bindings
option(BUILD_PYTHON "Build Python bindings for amsimpy" OFF)

if (BUILD_PYTHON)
  set(PYBIND11_FINDPYTHON ON)
  set(PYBIND11_DETAILED_ERROR_MESSAGES ON CACHE BOOL "Detailed error messages")
  find_package(pybind11 REQUIRED)

  pybind11_add_module(_core
      ${CMAKE_SOURCE_DIR}/python/bindings/python_bindings.cpp
      ${CMAKE_SOURCE_DIR}/python/bindings/context_wrapper.cpp)
  target_link_libraries(_core PRIVATE amsimcpp)
  target_include_directories(_core PRIVATE ${CMAKE_SOURCE_DIR}/include)
  target_compile_definitions(_core PRIVATE PYBIND11_DETAILED_ERROR_MESSAGES)
endif()

# Use GoogleTest

# include(FetchContent)
# set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
# set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)
# FetchContent_Declare(
#   googletest
#   GIT_REPOSITORY https://github.com/google/googletest.git
#   GIT_TAG v1.17.0
# )
# FetchContent_MakeAvailable(googletest)
#
# enable_testing()
#
# include(GoogleTest)

# Examples
option(BUILD_EXAMPLES "Build example applications" ON)

if (BUILD_EXAMPLES)
  add_executable(demo_simulation examples/demo_simulation.cc)
  target_link_libraries(demo_simulation PRIVATE amsimcpp)
endif()

# Install rules (so other projects can find+link amsimcpp)
include(GNUInstallDirs)
install(TARGETS amsimcpp
  EXPORT amsimcppTargets)
install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(EXPORT amsimcppTargets
  NAMESPACE amsimcpp::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/amsimcpp)
