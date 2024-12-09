cmake_minimum_required(VERSION 3.22)
project(Collinearity)

set(CMAKE_CXX_STANDARD 17)

include_directories(external/parlaylib/include)
include_directories(external/libvbyte)
include_directories(external/hopscotch-map/include)
include_directories(external/kseqpp/include)

add_subdirectory(external/streamvbyte)
add_subdirectory(external/kseqpp)

add_library(collinear STATIC
        prelude.h
        index_refs.cpp
        collinearity.h
        external/libvbyte/varintdecode.c
        external/libvbyte/vbyte.cc
        compressed_array.h
        parlay_utils.h
        cqueue.h
        cqutils.h)

add_executable(Collinearity main.cpp
        query.c)
add_executable(tests
        tests/testmain.cpp
        tests/test_compressed_array.cpp
        tests/test_cqueue.cpp
        tests/test_cqutils.cpp
        tests/tests.h)

target_compile_options(collinear PRIVATE -fpermissive -mavx)
target_compile_options(Collinearity PRIVATE -fpermissive)
target_link_libraries(Collinearity collinear z)

target_compile_options(tests PRIVATE -fpermissive)
target_link_libraries(tests collinear z)
