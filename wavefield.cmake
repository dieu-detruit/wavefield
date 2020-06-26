find_package(FFTW3 REQUIRED)

find_library (
    FFTW3_LIB
    NAMES fftw3 libfftw3 # what to look for
    HINTS ${FFTW3_LIBRARY_DIRS} # where to look
    NO_DEFAULT_PATH # do not search system default paths
)

find_library (
    FFTW3_THREAD_LIB
    NAMES fftw3_threads libfftw3_threads # what to look for
    HINTS ${FFTW3_LIBRARY_DIRS} # where to look
    NO_DEFAULT_PATH # do not search system default paths
)

find_package(TBB REQUIRED)

set(WAVEFIELD_INCLUDE_DIRS
        ${CMAKE_CURRENT_LIST_DIR}/include
        ${CMAKE_CURRENT_LIST_DIR}/external/grid/include
        ${FFTW3_INCLUDE_DIRS}
        )

set(WAVEFIELD_DEPENDENT_LIBS 
        ${FFTW3_LIB} ${FFTW3_THREAD_LIB} ${TBB_IMPORTED_TARGETS}
        )
