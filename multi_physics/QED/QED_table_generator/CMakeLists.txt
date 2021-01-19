
set(name table_generator)

add_executable(${name} "${name}.cpp")

target_link_libraries(${name} PRIVATE PXRMP_QED)

# OpenMP support
if(PXRMP_QED_OMP)
   find_package(OpenMP REQUIRED)
   target_link_libraries(${name} PRIVATE OpenMP::OpenMP_CXX)
   target_compile_definitions(${name} PRIVATE PXRMP_TABLE_GEN_HAS_OPENMP=1)
endif()

configure_file(table_inspector.ipynb table_inspector.ipynb COPYONLY)


# Move table generator and verification script in a tools subdirectory
set_target_properties(${name} PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/tools)

# Enforce C++14
target_compile_features(${name} PUBLIC cxx_std_14)
set_target_properties(${name} PROPERTIES CXX_EXTENSIONS OFF)

# Enable warnings
if(MSVC)
    target_compile_options(${name} PRIVATE /W4)
else()
    target_compile_options(${name} PRIVATE -Wall -Wextra -pedantic)
endif()
