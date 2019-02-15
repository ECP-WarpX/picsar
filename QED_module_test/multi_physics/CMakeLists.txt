cmake_minimum_required(VERSION 3.1)

set(CMAKE_CXX_EXTENSIONS OFF)
add_compile_options(-Wall -Wextra -Wpedantic)

add_library (multi_physics landau_lifshitz.cpp landau_lifshitz.h commons.h modified_bessel_functions.h modified_bessel_functions.cpp)

target_include_directories(multi_physics PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_features(multi_physics PUBLIC cxx_std_11)
