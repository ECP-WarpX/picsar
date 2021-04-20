
set(name kokkos_example)

add_executable(${name} "${name}.cpp")

target_link_libraries(${name} PRIVATE PXRMP_QED)


option(PXRMP_kokkos_internal                    "Download & build Kokkos" ON)

if(NOT PXRMP_KOKKOS_EXAMPLE)
    mark_as_advanced(PXRMP_kokkos_internal)
endif()

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

target_link_libraries (${name} PRIVATE Kokkos::kokkos)