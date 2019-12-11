cmake_minimum_required(VERSION 3.5)

set(CMAKE_CXX_EXTENSIONS OFF)

###Header only library called "pxrmp_qed"
add_library(pxrmp_qed INTERFACE)

###Finds Boost library (multithreaded version).
set(Boost_USE_MULTITHREADED      ON)
find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

target_include_directories(pxrmp_qed INTERFACE src/)
target_include_directories(pxrmp_qed INTERFACE src/containers)
target_include_directories(pxrmp_qed INTERFACE src/math)
target_include_directories(pxrmp_qed INTERFACE src/physics)
