add_library(QED INTERFACE)
add_library(${PROJECT_NAME}::QED ALIAS QED)

option(pxrmp_qed_test "Build PICSAR QED tests " ${IS_TOPLEVEL})

target_compile_features(QED INTERFACE cxx_std_14)

find_package(Boost 1.65.1)
if(Boost_FOUND)
    target_link_libraries(QED INTERFACE Boost::Boost)
    target_compile_definitions(QED INTERFACE PXRMP_HAS_BOOST=1)
endif()

target_include_directories(QED INTERFACE src/)
target_include_directories(QED INTERFACE src/containers)
target_include_directories(QED INTERFACE src/math)
target_include_directories(QED INTERFACE src/utils)
target_include_directories(QED INTERFACE src/physics)
target_include_directories(QED INTERFACE src/physics/schwinger)
target_include_directories(QED INTERFACE src/physics/breit_wheeler)
target_include_directories(QED INTERFACE src/physics/quantum_sync)
