cmake_minimum_required(VERSION 3.16)
project(amsimcpp LANGUAGES CXX VERSION 0.1.0)

# C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Library target
add_library(amsimcpp STATIC
  src/genome.cpp
  src/phenotype.cpp
  src/utils.cpp
)

# Public headers for users of the lib
target_include_directories(amsimcpp
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>)

# Reasonable warnings (optional)
if (CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
  target_compile_options(amsimcpp PRIVATE -Wall -Wextra -Wpedantic)
elseif (MSVC)
  target_compile_options(amsimcpp PRIVATE /W4)
endif()

# Example program (optional)
add_executable(demo examples/demo.cpp)
target_link_libraries(demo PRIVATE amsimcpp)

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
