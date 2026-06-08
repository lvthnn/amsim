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
else()
  set(BLA_STATIC ON)
  find_package(BLAS REQUIRED)
  find_package(LAPACK REQUIRED)
  target_link_libraries(amsimcpp INTERFACE BLAS::BLAS LAPACK::LAPACK)
endif()

# HDF5 support
set(HDF5_USE_STATIC_LIBRARIES ON)
find_package(HDF5 REQUIRED COMPONENTS C)
target_link_libraries(amsimcpp INTERFACE HDF5::HDF5)

# Fetch dependencies
include(FetchContent)

# Verbose downloads
set(FETCHCONTENT_QUIET False)

# Boost
set(BOOST_INCLUDE_LIBRARIES algorithm math)
set(BOOST_ENABLE_CMAKE ON)

FetchContent_Declare(
  Boost
  URL https://github.com/boostorg/boost/releases/download/boost-1.91.0-1/boost-1.91.0-1-cmake.tar.gz
  USES_TERMINAL_DOWNLOAD TRUE
  DOWNLOAD_EXTRACT_TIMESTAMP ON
  DOWNLOAD_NO_EXTRACT FALSE
  GIT_PROGRESS TRUE
)
FetchContent_MakeAvailable(Boost)

target_link_libraries(amsimcpp INTERFACE Boost::math Boost::algorithm)

# Eigen
FetchContent_Declare(
  Eigen
  GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
  GIT_TAG 5.0.1
  USES_TERMINAL_DOWNLOAD TRUE
  DOWNLOAD_EXTRACT_TIMESTAMP ON
  EXCLUDE_FROM_ALL
)
FetchContent_MakeAvailable(Eigen)

target_link_libraries(amsimcpp INTERFACE Eigen3::Eigen)

# HighFive
FetchContent_Declare(
  HighFive
  GIT_REPOSITORY https://github.com/highfive-devs/highfive.git
  GIT_TAG v3.3.0
  USES_TERMINAL_DOWNLOAD TRUE
  DOWNLOAD_EXTRACT_TIMESTAMP ON
  EXCLUDE_FROM_ALL
)
set(HIGHFIVE_USE_EIGEN ON CACHE BOOL "" FORCE)
set(HIGHFIVE_UNIT_TESTS OFF CACHE BOOL "" FORCE)
set(HIGHFIVE_BUILD_DOCS OFF CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(HighFive)

target_link_libraries(amsimcpp INTERFACE HighFive::HighFive)

# tomlplusplus
FetchContent_Declare(
    tomlplusplus
    GIT_REPOSITORY https://github.com/marzer/tomlplusplus.git
    GIT_TAG v3.4.0
    USES_TERMINAL_DOWNLOAD TRUE
    DOWNLOAD_EXTRACT_TIMESTAMP ON
)
FetchContent_MakeAvailable(tomlplusplus)

target_link_libraries(amsimcpp INTERFACE tomlplusplus::tomlplusplus)

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

target_compile_definitions(amsimcpp INTERFACE EIGEN_USE_BLAS EIGEN_USE_LAPACK)
target_compile_definitions(amsimcpp INTERFACE HIGHFIVE_USE_EIGEN)

# The CLI tool
option(BUILD_EXECUTABLE "Build the amsim CLI tool" ON)

if (BUILD_EXECUTABLE)
  add_executable(amsim src/amsim.cc)
  target_link_libraries(amsim PRIVATE amsimcpp)
  target_compile_definitions(amsim PRIVATE
    AMSIM_VERSION="${PROJECT_VERSION}"
    AMSIM_ARCH="${CMAKE_SYSTEM_PROCESSOR}"
    AMSIM_BUILD_TYPE="${CMAKE_BUILD_TYPE}"
    AMSIM_COMPILER="${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}"
  )
endif()
