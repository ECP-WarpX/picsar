cmake_minimum_required(VERSION 3.5)

set(CMAKE_CXX_EXTENSIONS OFF)

###Header only library called "pxrmp_qed"
add_library(pxrmp_qed INTERFACE)

###Finds Boost library (multithreaded version).
set(Boost_USE_MULTITHREADED      ON)
find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

###If the user specifies Kokkos directory, build with Kokkos support
if(WITH_KOKKOS)
  target_compile_definitions(pxrmp_qed INTERFACE PXRMP_BUILD_WITH_KOKKOS_SUPPORT)
  target_link_libraries(pxrmp_qed INTERFACE kokkos)
endif()

target_include_directories(pxrmp_qed INTERFACE src/)
