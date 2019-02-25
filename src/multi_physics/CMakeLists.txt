cmake_minimum_required (VERSION 3.1)

set(CMAKE_CXX_EXTENSIONS OFF)

#########THIS IS ONLY FOR DEBUG PURPOSES#########################################
add_compile_options(-Wall -Wextra -Wpedantic)
#################################################################################

project (testbed)

set(
    CMAKE_RUNTIME_OUTPUT_DIRECTORY
    ${CMAKE_BINARY_DIR}/bin
    )

add_subdirectory(QED)
add_subdirectory(QED_Tests)
