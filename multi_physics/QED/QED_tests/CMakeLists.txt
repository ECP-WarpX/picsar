# Check if unit_test package is available
find_package(Boost 1.66.0 REQUIRED COMPONENTS unit_test_framework)

set(TEST_NAMES
    picsar_algo
    picsar_array
    picsar_breit_wheeler_core
    picsar_breit_wheeler_tables
    picsar_breit_wheeler_tables_generator
    picsar_breit_wheeler_tabulated_functions
    picsar_chi_functions
    picsar_gamma_functions
    picsar_cmath_overload
    picsar_math_constants
    picsar_phys_constants
    picsar_quadrature
    picsar_quantum_sync_core
    picsar_quantum_sync_tables
    picsar_quantum_sync_tables_generator
    picsar_quantum_sync_tabulated_functions
    picsar_schwinger_engine_core
    picsar_serialization
    picsar_span
    picsar_spec_functions
    picsar_tables
    picsar_units
    picsar_vec_functions
    progress_bar
)

set(TEST_SRCS ${TEST_NAMES})
list(TRANSFORM TEST_SRCS PREPEND "test_")
list(TRANSFORM TEST_SRCS APPEND ".cpp")

# Run through each source
foreach(name ${TEST_NAMES})
    set(src "test_${name}.cpp")

    # Add compile target
    add_executable(${name} ${src})

    # link to Boost libraries AND your targets and dependencies
    target_link_libraries(${name} PRIVATE PXRMP_QED)
    target_link_libraries(${name} PRIVATE Boost::boost Boost::unit_test_framework)
    if(PXRMP_BOOST_TEST_DYN)
        target_compile_definitions(${name} PRIVATE BOOST_TEST_DYN_LINK)
    endif()

    # Move testing binaries in a bin subdirectory
    set_target_properties(${name} PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

    # Require C++14 or newer
    target_compile_features(${name} PUBLIC cxx_std_14)
    set_target_properties(${name} PROPERTIES CXX_EXTENSIONS OFF)

    # Enable warnings
    if(MSVC)
        target_compile_options(${name} PRIVATE /W4)
    else()
        target_compile_options(${name} PRIVATE -Wall -Wextra -pedantic)
    endif()


    # Finally add it to test execution
    message(STATUS "Test: ${name} (${src})")
    add_test(NAME ${name}
        COMMAND ${name}
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
endforeach()
