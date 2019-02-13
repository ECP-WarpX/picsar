cmake_minimum_required(VERSION 3.1)

set(CMAKE_CXX_EXTENSIONS OFF)
add_compile_options(-Wall -Wextra -Wpedantic)

add_library (picsar_multi_physics landau_lifshitz.cpp landau_lifshitz.h commons.h)

target_include_directories(picsar_multi_physics PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_features(testbed PUBLIC cxx_std_11)
