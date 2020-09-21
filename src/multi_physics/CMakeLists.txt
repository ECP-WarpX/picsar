cmake_minimum_required (VERSION 3.5)

set(CMAKE_CXX_EXTENSIONS OFF)

#Setup CMake to run tests
enable_testing()

project (qed)

add_subdirectory(QED)

add_definitions(-DBOOST_TEST_DYN_LINK)
add_subdirectory(QED_tests)
