cmake_minimum_required(VERSION 3.1)

set(CMAKE_CXX_EXTENSIONS OFF)

#########THIS IS ONLY FOR DEBUG PURPOSES#########################################
add_compile_options(-Wall -Wextra -Wpedantic)
#################################################################################


set(Boost_USE_MULTITHREADED      ON)
find_package(Boost COMPONENTS math_c99)
if(Boost_FOUND)
  add_definitions(-DHAS_BOOST_MATH)
endif()

add_library (multi_physics
            src/commons.h
            src/landau_lifshitz.cpp src/landau_lifshitz.h
            src/special_functions.h src/special_functions.cpp
            src/chi_calc_functions.h src/chi_calc_functions.cpp
            src/nonlin_breit_wheeler_engine.h src/nonlin_breit_wheeler_engine.cpp
            src/quadrature.h src/quadrature.cpp
            src/rng_wrapper.h src/rng_wrapper.cpp
            )

target_include_directories(multi_physics PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)
target_compile_features(multi_physics PUBLIC cxx_std_11)
