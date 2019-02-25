cmake_minimum_required(VERSION 3.0)

set(CMAKE_CXX_EXTENSIONS OFF)

#########THIS IS ONLY FOR DEBUG PURPOSES#########################################
add_compile_options(-Wall -Wextra -Wpedantic)
#################################################################################

project (testbed)
add_executable(testbed src/main.cpp src/species.h src/species.cpp src/photons.h
                        src/photons.cpp src/leptons.h  src/leptons.cpp  src/electrons.cpp
                        src/positrons.cpp src/electrons.h  src/positrons.h
                        )

target_compile_features(testbed PUBLIC cxx_std_11)

target_link_libraries(testbed qed)
