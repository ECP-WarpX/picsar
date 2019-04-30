cmake_minimum_required(VERSION 3.5)

#Check if unit_test package is available
find_package(Boost COMPONENTS unit_test_framework REQUIRED)
include_directories (${Boost_INCLUDE_DIRS})

file(GLOB DEMO_SRCS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cpp)

#Run through each source
foreach(demo_src ${DEMO_SRCS})
        #Extract the filename without an extension (NAME_WE)
        get_filename_component(demo_name ${demo_src} NAME_WE)

        #Add compile target
        add_executable(${demo_name} ${demo_src})

        #link to Boost libraries AND your targets and dependencies
        target_link_libraries(${demo_name} ${Boost_LIBRARIES} pxrmp_qed)

        #Move testing binaries in a demo_bin subdirectory
        set_target_properties(${test_name} PROPERTIES
            RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/demo_bin)

        #Enforce C++11
        set(CMAKE_CXX_EXTENSIONS OFF)
        if(NOT ${CMAKE_VERSION} VERSION_LESS "3.8.0")
          target_compile_features(${demo_name} PUBLIC cxx_std_11)
        else()
          set_target_properties(${demo_name} PROPERTIES CXX_STANDARD 11)
        endif()
endforeach(demo_src)
