#Check if unit_test package is available
find_package(Boost 1.65.1 COMPONENTS unit_test_framework REQUIRED)

file(GLOB TEST_SRCS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cpp)

#Run through each source
foreach(test_src ${TEST_SRCS})
        #Extract the filename without an extension (NAME_WE)
        get_filename_component(test_name ${test_src} NAME_WE)

        #Add compile target
        add_executable(${test_name} ${test_src})

        #link to Boost libraries AND your targets and dependencies
        target_link_libraries(${test_name} PRIVATE  ${Boost_LIBRARIES} pxrmp_qed)

        #Move testing binaries in a test_bin subdirectory
        set_target_properties(${test_name} PROPERTIES
            RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/test_bin)

        #Enforce C++14
        set(CMAKE_CXX_EXTENSIONS OFF)
        target_compile_features(${test_name} PUBLIC cxx_std_14)
        target_link_libraries(${test_name} PRIVATE Boost::unit_test_framework)
        target_link_libraries(${test_name} PRIVATE Boost::unit_test_framework)

        #Finally add it to test execution
        message(${test_name})
        message(${CMAKE_BINARY_DIR}/test_bin/${test_name})
        add_test(NAME ${test_name} COMMAND ${CMAKE_BINARY_DIR}/test_bin/${test_name})
endforeach(test_src)
