set(EXAMPLE_NAMES
    breit_wheeler
    quantum_synchrotron
)

# This macro downloads Kokkos from git if a PXRMP_kokkos_src is not provided
macro(find_kokkos)

    if(PXRMP_kokkos_src)
        message(STATUS "Compiling local Kokkos ...")
        message(STATUS "Kokkos source path: ${PXRMP_kokkos_src}")
    elseif(PXRMP_kokkos_internal)
        message(STATUS "Downloading Kokkos ...")
        message(STATUS "Kokkos repository: ${PXRMP_kokkos_repo} (${PXRMP_kokkos_tag})")
        include(FetchContent)
    else()
        message(FATAL_ERROR "Either PXRMP_kokkos_src or PXRMP_kokkos_internal must be defined!")
    endif()

    if(PXRMP_kokkos_src)
        list(APPEND CMAKE_MODULE_PATH "${PXRMP_kokkos_src}")
        add_subdirectory(${PXRMP_kokkos_src} _deps/localkokkos-build/)
    else()
        FetchContent_Declare(fetchedkokkos
        GIT_REPOSITORY ${PXRMP_kokkos_repo}
        GIT_TAG        ${PXRMP_kokkos_tag}
        BUILD_IN_SOURCE 0)
        FetchContent_GetProperties(fetchedkokkos)

        if(NOT fetchedkokkos_POPULATED)
            FetchContent_Populate(fetchedkokkos)
            list(APPEND CMAKE_MODULE_PATH "${fetchedkokkos_SOURCE_DIR}")
            add_subdirectory(${fetchedkokkos_SOURCE_DIR} ${fetchedkokkos_BINARY_DIR})
        endif()

        mark_as_advanced(FETCHCONTENT_BASE_DIR)
        mark_as_advanced(FETCHCONTENT_FULLY_DISCONNECTED)
        mark_as_advanced(FETCHCONTENT_QUIET)
        mark_as_advanced(FETCHCONTENT_SOURCE_DIR_FETCHEDKOKKOS)
        mark_as_advanced(FETCHCONTENT_UPDATES_DISCONNECTED)
        mark_as_advanced(FETCHCONTENT_UPDATES_DISCONNECTED_FETCHEDKOKKOS)

    endif()

endmacro()

# By disabling this the user can provide a path to a Kokkos source directory
option(PXRMP_kokkos_internal                    "Download & build Kokkos" ON)

if(NOT PXRMP_KOKKOS_EXAMPLE)
    mark_as_advanced(PXRMP_kokkos_internal)
endif()

# local source-tree
set(PXRMP_kokkos_src ""
    CACHE PATH
    "Local path to Kokkos source directory (preferred if set)")

# Git fetcher
set(PXRMP_kokkos_repo "https://github.com/kokkos/kokkos.git"
    CACHE STRING
    "Repository URI to pull and build AMReX from if(PXRMP_kokkos_internal)")
set(PXRMP_kokkos_tag "3.3.01"
    CACHE STRING
    "git tag for PXRMP_kokkos_repo if(PXRMP_kokkos_internal)")

find_kokkos()


# Run through each source
foreach(name ${EXAMPLE_NAMES})
    set(src "kokkos_example_${name}.cpp")

    # Add compile target
    add_executable(${name} ${src})

    target_link_libraries(${name} PRIVATE PXRMP_QED)
    target_link_libraries (${name} PRIVATE Kokkos::kokkos)

    # Move examples into /bin/kokkos_examples
    set_target_properties(${name} PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/kokkos_examples)

    # Require C++14 or newer
    target_compile_features(${name} PUBLIC cxx_std_14)
    set_target_properties(${name} PROPERTIES CXX_EXTENSIONS OFF)

    message(STATUS "Example: ${name} (${src})")

endforeach()
