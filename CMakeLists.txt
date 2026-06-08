cmake_minimum_required(VERSION 3.10)
project(anewdsc C)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Wextra -g")

# Find GMP, MPFR, MPFI
find_library(GMP_LIBRARY gmp REQUIRED)
find_library(MPFR_LIBRARY mpfr REQUIRED)
find_library(MPFI_LIBRARY mpfi REQUIRED)

include_directories(include)

add_library(anewdsc STATIC
    src/anewdsc_poly.c
    src/anewdsc_roots.c
)

target_link_libraries(anewdsc
    ${GMP_LIBRARY}
    ${MPFR_LIBRARY}
    ${MPFI_LIBRARY}
)
target_include_directories(anewdsc PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

enable_testing()
add_executable(test_poly tests/test_poly.c)
target_link_libraries(test_poly anewdsc)

# Examples
add_executable(random_roots examples/random_roots.c)
target_link_libraries(random_roots anewdsc)
add_test(NAME test_poly COMMAND test_poly)
