cmake_minimum_required(VERSION 3.10)
project(SkewOptimizer)

find_package(PkgConfig REQUIRED)
pkg_check_modules(GMP REQUIRED IMPORTED_TARGET gmp)

add_executable(skewopt skewopt.cpp msieve_poly.cpp)
target_link_libraries(skewopt PRIVATE PkgConfig::GMP)