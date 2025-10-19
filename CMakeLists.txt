cmake_minimum_required(VERSION 3.20)
project(amsimcpp LANGUAGES CXX VERSION 0.1.0)

# C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_BUILD_TYPE Release) # for optimisation and speedup

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
  src/stats.cc
)

# Public headers for users of the lib
target_include_directories(amsimcpp
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>)

if (APPLE)
  find_library(ACCELERATE_FRAMEWORK Accelerate REQUIRED)
  target_link_libraries(amsimcpp PUBLIC ${ACCELERATE_FRAMEWORK})
  target_compile_definitions(amsimcpp PUBLIC ACCELERATE_NEW_LAPACK)
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

# Reasonable warnings (optional)
if (CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
  target_compile_options(amsimcpp PRIVATE -Wall -Wextra -Wpedantic)
elseif (MSVC)
  target_compile_options(amsimcpp PRIVATE /W4)
endif()

# Example program (optional)
add_executable(demo_simulation examples/demo_simulation.cc)
add_executable(demo_stats      examples/demo_stats.cc)

target_link_libraries(demo_simulation PRIVATE amsimcpp)
target_link_libraries(demo_stats PRIVATE amsimcpp)

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
