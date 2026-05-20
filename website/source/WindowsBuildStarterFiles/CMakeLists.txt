cmake_minimum_required(VERSION 3.10)
project(cownoise)

find_package(PkgConfig REQUIRED)
pkg_check_modules(GMP REQUIRED IMPORTED_TARGET gmp)

add_executable(cyclo cyclo.cpp)
target_link_libraries(cyclo PRIVATE PkgConfig::GMP)

add_executable(cyclop cyclop.cpp)
target_link_libraries(cyclop PRIVATE PkgConfig::GMP)
