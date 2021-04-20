set(name pxr_qed)

set(_PY_DEV_MODULE Development.Module)
if(CMAKE_VERSION VERSION_LESS 3.18.0)
    # over-specification needed for CMake<3.18
    #   https://pybind11.readthedocs.io/en/latest/compiling.html#findpython-mode
    #   https://cmake.org/cmake/help/v3.18/module/FindPython.html
    set(_PY_DEV_MODULE Development)
endif()
find_package(Python COMPONENTS Interpreter ${_PY_DEV_MODULE} REQUIRED)
find_package(pybind11 2.6.0 CONFIG REQUIRED)

pybind11_add_module(${name} "${name}.cpp")

target_link_libraries(${name} PRIVATE PXRMP_QED)

option(
    PXRQEDPY_DOUBLE_PRECISION
    "Compile pxr_qed python module in double precision" ON)

if(PXRQEDPY_DOUBLE_PRECISION)
    target_compile_definitions(${name} PRIVATE PXRQEDPY_DOUBLE_PRECISION=1)
endif()

set(PXRQEDPY_UNITS_VALUES SI NORM_OMEGA NORM_LAMBDA HEAVISIDE_LORENTZ)
set(PXRQEDPY_UNITS SI CACHE STRING "Unit system for the pxr_qed python module (SI/NORM_OMEGA/NORM_LAMBDA/HEAVISIDE_LORENTZ)")
set_property(CACHE PXRQEDPY_UNITS PROPERTY STRINGS ${PXRQEDPY_UNITS_VALUES})

if(NOT PXRQEDPY_UNITS IN_LIST PXRQEDPY_UNITS_VALUES)
    message(FATAL_ERROR "PXRQEDPY_UNITS (${PXRQEDPY_UNITS}) must be one of ${PXRQEDPY_UNITS_VALUES}")
endif()

if(PXRQEDPY_UNITS STREQUAL "SI")
    target_compile_definitions(${name} PRIVATE PXRQEDPY_SI=1)
elseif(PXRQEDPY_UNITS STREQUAL "NORM_OMEGA")
    target_compile_definitions(${name} PRIVATE PXRQEDPY_NORM_OMEGA=1)
elseif(PXRQEDPY_UNITS STREQUAL "NORM_LAMBDA")
    target_compile_definitions(${name} PRIVATE PXRQEDPY_NORM_LAMBDA=1)
elseif(PXRQEDPY_UNITS STREQUAL "HEAVISIDE_LORENTZ")
    target_compile_definitions(${name} PRIVATE PXRQEDPY_HEAVISIDE_LORENTZ=1)
else()
    message(FATAL_ERROR "PXRQEDPY_UNITS is incorrectly defined")
endif()

# OpenMP support
if(PXRMP_QED_OMP)
   find_package(OpenMP REQUIRED)
   target_link_libraries(${name} PRIVATE OpenMP::OpenMP_CXX)
   target_compile_definitions(${name} PRIVATE PXQEDPY_HAS_OPENMP=1)
endif()

# Move module in "python_bindings" subdirectory
set_target_properties(${name} PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/python_bindings)

configure_file(demo_python_bindings.ipynb
    ${CMAKE_BINARY_DIR}/python_bindings/demo_python_bindings.ipynb COPYONLY)

# Require C++14 or newer
target_compile_features(${name} PUBLIC cxx_std_14)
set_target_properties(${name} PROPERTIES CXX_EXTENSIONS OFF)

# Enable warnings
if(MSVC)
    target_compile_options(${name} PRIVATE /W4)
else()
    target_compile_options(${name} PRIVATE -Wall -Wextra -pedantic)
endif()
