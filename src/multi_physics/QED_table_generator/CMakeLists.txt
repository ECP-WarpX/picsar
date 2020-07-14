cmake_minimum_required (VERSION 3.5)

set(CMAKE_CXX_EXTENSIONS OFF)

project (qed_table_generator)

set(
    CMAKE_RUNTIME_OUTPUT_DIRECTORY
    ${CMAKE_BINARY_DIR}/bin
    )

include_directories(../QED)

add_executable(qed_table_generator table_generator.cpp)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(qed_table_generator PUBLIC OpenMP::OpenMP_CXX)
endif()

#Optimize build for machine architecture 
option(OPTIMIZE_FOR_NATIVE "If supported, build with -march=native" OFF)
if(OPTIMIZE_FOR_NATIVE)
    include(CheckCXXCompilerFlag)
    check_cxx_compiler_flag("-march=native" COMPILER_SUPPORTS_NATIVE)
    if(COMPILER_SUPPORTS_NATIVE)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
    endif()
endif()
#_______________________________________

