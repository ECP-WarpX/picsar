cmake_minimum_required (VERSION 3.5)

set(CMAKE_CXX_EXTENSIONS OFF)

#Setup CMake to run tests
enable_testing()

project (qed)

###If the user specifies Kokkos directory, build with Kokkos support
option(WITH_KOKKOS "Build with Kokkos. User should provide KOKKOS_DIRECTORY" OFF)
if(WITH_KOKKOS)
  add_subdirectory(${KOKKOS_DIRECTORY}  ${CMAKE_BINARY_DIR}/kokkos)
  include_directories(${Kokkos_INCLUDE_DIRS_RET})
  message("The library will be built with Kokkos support")
endif()

add_subdirectory(QED)
add_subdirectory(QED_tests)
