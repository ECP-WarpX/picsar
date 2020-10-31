find_package(Boost)
if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
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
