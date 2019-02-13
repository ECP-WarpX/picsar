cmake_minimum_required(VERSION 3.1)

set(CMAKE_CXX_STANDARD 11)

set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pedantic ")

project (testbed)
add_executable(testbed src/main.cpp src/commons.h src/species.h src/species.cpp src/photons.h
                        src/photons.cpp src/leptons.h  src/leptons.cpp  src/electrons.cpp
                        src/positrons.cpp src/electrons.h  src/positrons.h)
