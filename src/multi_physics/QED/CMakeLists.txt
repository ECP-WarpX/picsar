add_library(QED INTERFACE)
add_library(${PROJECT_NAME}::QED ALIAS QED)

target_compile_features(QED INTERFACE cxx_std_14)

if(PXRMP_QED_OMP)
   find_package(OpenMP REQUIRED)
   target_link_libraries(QED INTERFACE OpenMP::OpenMP_CXX)
   target_compile_definitions(QED INTERFACE PXRMP_HAS_OPENMP=1)
endif()

find_package(Boost 1.66.0)
if(Boost_FOUND)
    target_link_libraries(QED INTERFACE Boost::boost)
    target_compile_definitions(QED INTERFACE PXRMP_HAS_BOOST=1)
endif()

if(PXRMP_QED_TABLEGEN)
    if(NOT Boost_FOUND)
         message( FATAL_ERROR "Table generation requires Boost (>= 1.66)!" )
    endif()
endif()

target_include_directories(QED INTERFACE src/)
target_include_directories(QED INTERFACE src/containers)
target_include_directories(QED INTERFACE src/math)
target_include_directories(QED INTERFACE src/utils)
target_include_directories(QED INTERFACE src/physics)
target_include_directories(QED INTERFACE src/physics/schwinger)
target_include_directories(QED INTERFACE src/physics/breit_wheeler)
target_include_directories(QED INTERFACE src/physics/quantum_sync)

# Setup CMake to run tests
if(PXRMP_QED_TEST)
  add_subdirectory(tests)
endif()

