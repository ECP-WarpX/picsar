cmake_minimum_required(VERSION 3.12)

project(
  PXRMP_QED
  DESCRIPTION
    "A header only C++ library (part of the PICSAR library) providing QED modules for Particle-In-Cell codes."
  HOMEPAGE_URL "https://github.com/ECP-WarpX/picsar")

if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
  set(IS_TOPLEVEL TRUE)
else()
  set(IS_TOPLEVEL FALSE)
endif()


########
# Options
########
option(PXRMP_QED_OMP              "Enable OpenMP support" ON)
option(PXRMP_QED_TABLEGEN         "Enable table generation (needs Boost)" ON)
option(PXRMP_QED_TEST             "Build PICSAR QED tests " ${IS_TOPLEVEL})
option(PXRMP_QED_TOOLS            "Build PICSAR QED tools " ${IS_TOPLEVEL})
option(PXRMP_BOOST_TEST_DYN_LINK  "Link against the Boost Unit Test Framework shared library" ON)
option(PXRMP_DPCPP_FIX            "Use cl::sycl::floor cl::sycl::floorf on device" OFF)


########
# Library include directory
########
add_library(PXRMP_QED INTERFACE)
add_library(${PROJECT_NAME}::PXRMP_QED ALIAS PXRMP_QED)

target_compile_features(PXRMP_QED INTERFACE cxx_std_14)

if(PXRMP_QED_OMP)
   find_package(OpenMP REQUIRED)
   target_link_libraries(PXRMP_QED INTERFACE OpenMP::OpenMP_CXX)
   target_compile_definitions(PXRMP_QED INTERFACE PXRMP_HAS_OPENMP=1)
endif()

find_package(Boost 1.66.0)
if(Boost_FOUND)
    target_link_libraries(PXRMP_QED INTERFACE Boost::boost)
    target_compile_definitions(PXRMP_QED INTERFACE PXRMP_HAS_BOOST=1)
endif()

if(PXRMP_QED_TABLEGEN)
    if(NOT Boost_FOUND)
         message( FATAL_ERROR "Table generation requires Boost (>= 1.66)!" )
    endif()
endif()

if(PXRMP_DPCPP_FIX)
    target_compile_definitions(PXRMP_QED INTERFACE PXRMP_DPCPP_FIX=1)
endif()

target_include_directories(
    ${PROJECT_NAME}
    INTERFACE
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include/>
    $<INSTALL_INTERFACE:include>
)


########
# Setup CMake to build tests
########
if(PXRMP_QED_TEST)
    if(NOT PXRMP_QED_TABLEGEN)
        message( FATAL_ERROR "Tests require to enable table generation!" )
    endif()
    enable_testing()
    add_subdirectory(QED_tests)
endif()


########
# Setup CMake to build tools
########
if(PXRMP_QED_TOOLS)
    if(NOT PXRMP_QED_TABLEGEN)
         message( FATAL_ERROR "Tools require to enable table generation!" )
    endif()
    add_subdirectory(QED_table_generator)
endif()
