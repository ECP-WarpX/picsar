cmake_minimum_required(VERSION 3.12)

project(
  pxrmp_qed
  DESCRIPTION
    "A header only c++ library (part of the PICSAR library) providing QED modules for Particle-In-Cell codes."
  HOMEPAGE_URL "https://github.com/ECP-WarpX/picsar"
  LANGUAGES CXX)
  
add_library(${PROJECT_NAME} INTERFACE)
add_library(${PROJECT_NAME}::${PROJECT_NAME} ALIAS ${PROJECT_NAME})

if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
  set(IS_TOPLEVEL TRUE)
else()
  set(IS_TOPLEVEL FALSE)
endif()

option(pxrmp_qed_test "Build PICSAR QED tests " ${IS_TOPLEVEL})

target_compile_features(${PROJECT_NAME} INTERFACE cxx_std_14)
add_subdirectory(QED)

#Setup CMake to run tests
if(pxrmp_qed_test)
  enable_testing()
  add_subdirectory(QED_tests)
endif()





