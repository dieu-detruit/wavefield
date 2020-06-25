cmake_minimum_required(VERSION 3.0)

set(CMAKE_CXX_COMPILER "g++-10" CACHE STRING "g++ compiler" FORCE)
set(CMAKE_C_COMPILER "gcc-10" CACHE STRING "gcc compiler" FORCE)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

project(Wavefield LANGUAGES C CXX)

include(wavefield.cmake)

# add test executables
file(GLOB TEST_SOURCES testsuite/*.cpp)
foreach(TEST_SOURCE ${TEST_SOURCES})
    get_filename_component(TEST_NAME ${TEST_SOURCE} NAME_WE)

    add_executable(${TEST_NAME} ${TEST_SOURCE})
    target_include_directories(${TEST_NAME} PRIVATE ${WAVEFIELD_INCLUDE_DIRS})
    target_link_libraries(${TEST_NAME} PRIVATE ${WAVEFIELD_DEPENDENT_LIBS})
    target_compile_options(${TEST_NAME} PRIVATE
            -g -O3)
endforeach()
