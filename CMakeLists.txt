cmake_minimum_required(VERSION 3.20)
project(amsimcpp LANGUAGES CXX VERSION 0.1.0)

# C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Add debugging options depending on CMAKE_BUILD_TYPE
if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    message(STATUS "Enabling sanitizers")
    add_compile_options(-fsanitize=address,undefined -fno-omit-frame-pointer)
    add_link_options(-fsanitize=address,undefined -fno-omit-frame-pointer)
endif()

# Library targets
add_library(amsimcpp STATIC
  src/genome.cpp
  src/phenotype.cpp
  src/utils.cpp
  src/mating.cpp
  src/phenobuf.cpp
  src/phenoarch.cpp
  src/haplobuf.cpp
  src/logger.cpp
  src/logger_timer.cpp
  src/metric.cpp
  src/metricspec.cpp
  src/simulation.cpp
  src/simulation_config.cpp
  src/simulation_results.cpp
  src/stats.cpp
)

# Public headers for users of the lib
target_include_directories(amsimcpp
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>)

# Resolve OS-dependent BLAS / LAPACK imports
if (APPLE)
  find_library(ACCELERATE_FRAMEWORK Accelerate REQUIRED)
  target_link_libraries(amsimcpp PUBLIC ${ACCELERATE_FRAMEWORK})
  target_compile_definitions(amsimcpp PUBLIC ACCELERATE_NEW_LAPACK)
  target_compile_options(amsimcpp PRIVATE -Wno-deprecated-declarations)
else()
  find_package(BLAS REQUIRED)
  find_package(LAPACK REQUIRED)

  if(TARGET BLAS::BLAS)
    target_link_libraries(amsimcpp PUBLIC BLAS::BLAS)
  else()
    target_link_libraries(amsimcpp PUBLIC ${BLAS_LIBRARIES})
  endif()

  if(TARGET LAPACK::LAPACK)
    target_link_libraries(amsimcpp PUBLIC LAPACK::LAPACK)
  else()
    target_link_libraries(amsimcpp PUBLIC ${LAPACK_LIBRARIES})
  endif()
endif()

# Python bindings
option(BUILD_PYTHON "Build Python bindings for amsimpy" OFF)

if (BUILD_PYTHON)
  set(PYBIND11_FINDPYTHON ON)
  find_package(pybind11 REQUIRED)

  pybind11_add_module(_core python/bindings/python_bindings.cpp)
  target_link_libraries(_core PRIVATE amsimcpp)
  target_include_directories(_core PRIVATE include)
endif()

# Reasonable warnings
if (CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
  target_compile_options(amsimcpp PRIVATE -Wall -Wextra -Wpedantic)
elseif (MSVC)
  target_compile_options(amsimcpp PRIVATE /W4)
endif()

# Use GoogleTest
include(FetchContent)
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)
FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG v1.17.0
)
FetchContent_MakeAvailable(googletest)

enable_testing()

# Example programs
add_executable(demo_simulation examples/demo_simulation.cpp)
add_executable(demo_manual examples/demo_manual.cpp)
add_executable(demo_results examples/demo_results.cpp)
add_executable(demo_logger examples/demo_logger.cpp)

target_link_libraries(demo_simulation PRIVATE amsimcpp)
target_link_libraries(demo_manual PRIVATE amsimcpp)
target_link_libraries(demo_results PRIVATE amsimcpp)
target_link_libraries(demo_logger PRIVATE amsimcpp)

# Link unit tests
add_executable(stats_test tests/stats_test.cpp)

target_link_libraries(stats_test PRIVATE GTest::gtest_main amsimcpp)

include(GoogleTest)
gtest_discover_tests(stats_test)


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
