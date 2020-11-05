cmake_minimum_required(VERSION 3.12)

project(
  PXRMP
  DESCRIPTION
    "A header only c++ library (part of the PICSAR library) providing QED modules for Particle-In-Cell codes."
  HOMEPAGE_URL "https://github.com/ECP-WarpX/picsar")

if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
  set(IS_TOPLEVEL TRUE)
else()
  set(IS_TOPLEVEL FALSE)
endif()

enable_testing()

option(PXRMP_QED_OMP      "Enable OpenMP support" ON)
option(PXRMP_QED_TABLEGEN "Enable table generation (needs Boost)" ON)
option(PXRMP_QED_TEST     "Build PICSAR QED tests " ${IS_TOPLEVEL})
option(PXRMP_QED_TOOLS    "Build PICSAR QED tools " ${IS_TOPLEVEL})

add_subdirectory(QED)
