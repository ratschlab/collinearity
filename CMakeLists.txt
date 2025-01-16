cmake_minimum_required(VERSION 3.22)
project(Collinearity LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)

include_directories(external/parlaylib/include)
include_directories(external/libvbyte)
include_directories(external/hopscotch-map/include)
include_directories(external/kseqpp/include)
include_directories(external/parlayhash/include)
include_directories(external/emhash)

add_subdirectory(external/streamvbyte)
add_subdirectory(external/kseqpp)
add_subdirectory(external/pybind11)

set(Python3_EXECUTABLE /home/sayan/.virtualenvs/SequenceAnalysis/bin/python3)
find_package(Python3 COMPONENTS Interpreter Development NumPy REQUIRED)
message(STATUS "Python3_EXECUTABLE: ${Python3_EXECUTABLE}")
message(STATUS "Python3_INCLUDE_DIRS: ${Python3_INCLUDE_DIRS}")
message(STATUS "Python3_LIBRARIES: ${Python3_LIBRARIES}")

set(PYBIND11_FINDPYTHON ON)
pybind11_add_module(pycollinearity pybindings.cpp)

add_library(collinear STATIC
        prelude.h
        external/libvbyte/varintdecode.c
        external/libvbyte/vbyte.cc
        compressed_array.h
        parlay_utils.h
        cqueue.h
        cqutils.h)

add_executable(Collinearity main.cpp
        index_refs.cpp
        query.cpp
        collinearity.h
        index.h
        xxhash.h
)
add_executable(tests
        tests/testmain.cpp
        tests/test_compressed_array.cpp
        tests/test_cqueue.cpp
        tests/test_cqutils.cpp
        tests/tests.h
)

target_compile_options(collinear PRIVATE -fpermissive -mavx)
target_compile_options(Collinearity PRIVATE -fpermissive)
target_link_libraries(Collinearity collinear z)

target_compile_options(tests PRIVATE -fpermissive)
target_link_libraries(tests collinear z)
