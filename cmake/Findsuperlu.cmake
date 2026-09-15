# Findsuperlu.cmake
#
# Locates the SuperLU sequential sparse direct solver.
#
# Search hints:
#   SUPERLU_ROOT / SUPERLU_DIR (variable or environment) - install prefix
#
# Defines:
#   SUPERLU_FOUND         - TRUE if SuperLU was located
#   SUPERLU_INCLUDE_DIRS  - include directories for SuperLU headers
#   SUPERLU_LIBRARIES     - libraries to link against
#   SUPERLU_LIBRARY_DIRS  - directory containing the SuperLU library
#
# Imported target:
#   superlu::superlu

set(SUPERLU_FOUND FALSE)

set(SUPERLU_INCLUDE_DIRS "")
set(SUPERLU_LIBRARY_DIRS "")
set(SUPERLU_LIBRARIES "")

find_path(SUPERLU_INCLUDE_DIR
    NAMES slu_ddefs.h supermatrix.h slu_util.h
    HINTS
        ${SUPERLU_ROOT}
        ${SUPERLU_DIR}
        ENV SUPERLU_ROOT
        ENV SUPERLU_DIR
        ${CMAKE_INSTALL_PREFIX}
    PATH_SUFFIXES include include/superlu superlu
)

find_library(SUPERLU_LIBRARY
    NAMES superlu superlu_5 superlu_6
    HINTS
        ${SUPERLU_ROOT}
        ${SUPERLU_DIR}
        ENV SUPERLU_ROOT
        ENV SUPERLU_DIR
        ${CMAKE_INSTALL_PREFIX}
    PATHS "${SUPERLU_INCLUDE_DIR}/../lib" "${SUPERLU_INCLUDE_DIR}/../../lib"
    PATH_SUFFIXES lib lib64
)

if(SUPERLU_INCLUDE_DIR AND SUPERLU_LIBRARY)
    set(SUPERLU_FOUND TRUE)
    set(SUPERLU_INCLUDE_DIRS ${SUPERLU_INCLUDE_DIR})
    set(SUPERLU_LIBRARIES ${SUPERLU_LIBRARY})
    get_filename_component(SUPERLU_LIBRARY_DIRS "${SUPERLU_LIBRARY}" DIRECTORY)

    # SuperLU calls into BLAS; a static libsuperlu.a will not resolve those
    # symbols on its own, so pull in BLAS as a dependency of the target.
    find_package(BLAS)
    if(BLAS_FOUND)
        if(TARGET BLAS::BLAS)
            list(APPEND SUPERLU_LIBRARIES BLAS::BLAS)
        else()
            list(APPEND SUPERLU_LIBRARIES ${BLAS_LIBRARIES})
        endif()
    endif()

    if(NOT TARGET superlu::superlu)
        add_library(superlu::superlu INTERFACE IMPORTED)
        set_target_properties(superlu::superlu PROPERTIES
            INTERFACE_LINK_LIBRARIES "${SUPERLU_LIBRARY}"
            INTERFACE_INCLUDE_DIRECTORIES "${SUPERLU_INCLUDE_DIRS}"
        )
    endif()
endif()

if(SUPERLU_FOUND)
    message(STATUS "Found SuperLU: ${SUPERLU_LIBRARIES}")
    message(STATUS "SuperLU include dir: ${SUPERLU_INCLUDE_DIRS}")
else()
    message(WARNING "SuperLU not found. Set SUPERLU_ROOT to the SuperLU install prefix.")
endif()

mark_as_advanced(SUPERLU_INCLUDE_DIR SUPERLU_LIBRARY)

set(SUPERLU_INCLUDE_DIRS ${SUPERLU_INCLUDE_DIRS} CACHE PATH "SuperLU include directories")
set(SUPERLU_LIBRARIES ${SUPERLU_LIBRARIES} CACHE PATH "SuperLU libraries")
set(SUPERLU_LIBRARY_DIRS ${SUPERLU_LIBRARY_DIRS} CACHE PATH "SuperLU library directories")
set(SUPERLU_FOUND ${SUPERLU_FOUND} CACHE BOOL "Is SuperLU found?")
