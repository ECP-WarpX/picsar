cmake_minimum_required (VERSION 3.5)

set(CMAKE_CXX_EXTENSIONS OFF)

project (testbed)

set(
    CMAKE_RUNTIME_OUTPUT_DIRECTORY
    ${CMAKE_BINARY_DIR}/bin_test_qed
    )

include_directories(../QED)

add_executable(testbed src/main.cpp src/species.h src/species.cpp src/photons.h
                        src/photons.cpp src/leptons.h  src/leptons.cpp  src/electrons.cpp
                        src/positrons.cpp src/electrons.h  src/positrons.h
                        )
