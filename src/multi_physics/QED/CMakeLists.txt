find_package(Boost 1.65.1 COMPONENTS math)
if(Boost_FOUND)
    target_link_libraries(pxrmp_qed INTERFACE Boost::Boost Boost::math)
    # TODO: define a preprocessor variable to let the library know if Boost is available or not
endif()

target_include_directories(${PROJECT_NAME} INTERFACE src/)
target_include_directories(${PROJECT_NAME} INTERFACE src/containers)
target_include_directories(${PROJECT_NAME} INTERFACE src/math)
target_include_directories(${PROJECT_NAME} INTERFACE src/utils)
target_include_directories(${PROJECT_NAME} INTERFACE src/physics)
target_include_directories(${PROJECT_NAME} INTERFACE src/physics/schwinger)
target_include_directories(${PROJECT_NAME} INTERFACE src/physics/breit_wheeler)
target_include_directories(${PROJECT_NAME} INTERFACE src/physics/quantum_sync)
