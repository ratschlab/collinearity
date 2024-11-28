cmake_minimum_required(VERSION 3.22)
project(Collinearity)

set(CMAKE_CXX_STANDARD 17)

include_directories(external/parlaylib/include)
include_directories(external/libvbyte)
include_directories(external/hopscotch-map/include)
include_directories(external/kseqpp/include)

add_subdirectory(external/streamvbyte)
add_subdirectory(external/kseqpp)

add_executable(Collinearity main.cpp
        prelude.h
        index_refs.cpp
        collinearity.h
        external/libvbyte/varintdecode.c
        external/libvbyte/vbyte.cc
        compressed_array.h
        parlay_utils.h
)

target_compile_options(Collinearity PRIVATE -fpermissive -mavx)

target_link_libraries(Collinearity z)
