cmake_minimum_required(VERSION 3.10)
project(phi)

find_package(PkgConfig REQUIRED)
pkg_check_modules(GMP REQUIRED IMPORTED_TARGET gmp)

add_executable(phi phi.c)
target_link_libraries(phi PRIVATE PkgConfig::GMP)
