# Check if unit_test package is available
find_package(Boost 1.66.0 REQUIRED COMPONENTS unit_test_framework)

set(TEST_NAMES
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
