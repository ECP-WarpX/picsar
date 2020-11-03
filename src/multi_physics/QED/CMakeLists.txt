add_library(QED INTERFACE)
add_library(${PROJECT_NAME}::QED ALIAS QED)

option(PXRMP_QED_OMP      "Enable OpenMP support" ON) # TODO: implement definitions
option(PXRMP_QED_TABLEGEN "Enable table generation (needs Boost)" ON)
option(PXRMP_QED_TEST     "Build PICSAR QED tests " ${IS_TOPLEVEL})

target_compile_features(QED INTERFACE cxx_std_14)

if(PXRMP_QED_OMP)
   find_package(OpenMP REQUIRED)
   target_link_libraries(QED INTERFACE OpenMP::OpenMP_CXX)
endif()

if(PXRMP_QED_TABLEGEN)
    find_package(Boost 1.66.0 REQUIRED)
    target_link_libraries(QED INTERFACE Boost::boost)
    target_compile_definitions(QED INTERFACE PXRMP_HAS_BOOST=1)
    target_compile_definitions(QED INTERFACE PXRMP_HAS_OPENMP=1) # FIXME: in-source definition guards
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

