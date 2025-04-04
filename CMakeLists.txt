cmake_minimum_required(VERSION 3.22)
project(Collinearity LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)

if ( NOT CMAKE_BUILD_TYPE )
    message(STATUS "Build Type not set, defaulting to Release..." )
    set( CMAKE_BUILD_TYPE Release )
else ()
    message(STATUS "Build type is ${CMAKE_BUILD_TYPE}" )
endif()

include_directories(external/parlaylib/include)
include_directories(external/kseqpp/include)
include_directories(external/emhash)
include_directories(external/argparse/include)
include_directories(external/libvbyte)

add_subdirectory(external/kseqpp)
add_subdirectory(external/pybind11)
add_subdirectory(external/slow5lib)

#set(Python3_EXECUTABLE /home/sayan/.virtualenvs/pycollinearity/bin/python3)
find_package(Python3 COMPONENTS Interpreter Development REQUIRED)
message(STATUS "Python3_EXECUTABLE: ${Python3_EXECUTABLE}")
message(STATUS "Python3_INCLUDE_DIRS: ${Python3_INCLUDE_DIRS}")
message(STATUS "Python3_LIBRARIES: ${Python3_LIBRARIES}")

set(PYBIND11_FINDPYTHON ON)
pybind11_add_module(_core
        src/index.cpp
        src/index_refs.cpp
        src/query.cpp
        src/pybindings.cpp
        external/libvbyte/varintdecode.c
        external/libvbyte/vbyte.cc
)

add_executable(Collinearity src/main.cpp
        src/index_refs.cpp
        src/query.cpp
        src/collinearity.h
        src/index.h
        src/rawsignals.h
        src/tstatsegmentation.cpp
        src/rawreferences.cpp
        src/rawsegmentation.cpp
        src/rawsignals.cpp
        src/config.h
        src/utils.h
        src/index.cpp
        src/mempool.h
        external/libvbyte/varintdecode.c
        external/libvbyte/vbyte.cc
)
#add_executable(test
#        scratch/testmain.cpp
#        scratch/tests.h
#)

#target_compile_options(collinear PRIVATE -fpermissive -mavx)
target_compile_options(Collinearity PRIVATE -fpermissive -mavx)
target_link_libraries(Collinearity slow5 z)

#target_compile_options(test PRIVATE -fpermissive)
#target_link_libraries(test z)

target_link_libraries(_core PUBLIC -mavx slow5 z pybind11::headers)
install(TARGETS _core DESTINATION pycollinearity)
