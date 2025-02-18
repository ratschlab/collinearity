cmake_minimum_required(VERSION 3.22)
project(Collinearity LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)

include_directories(external/parlaylib/include)
include_directories(external/kseqpp/include)
include_directories(external/emhash)
include_directories(external/argparse/include)

add_subdirectory(external/kseqpp)
add_subdirectory(external/pybind11)
add_subdirectory(external/slow5lib)

#set(Python3_EXECUTABLE /home/sayan/.virtualenvs/SequenceAnalysis/bin/python3)
find_package(Python3 COMPONENTS Interpreter Development REQUIRED)
message(STATUS "Python3_EXECUTABLE: ${Python3_EXECUTABLE}")
message(STATUS "Python3_INCLUDE_DIRS: ${Python3_INCLUDE_DIRS}")
message(STATUS "Python3_LIBRARIES: ${Python3_LIBRARIES}")

set(PYBIND11_FINDPYTHON ON)
pybind11_add_module(_core
        src/index.cpp
        src/index_refs.cpp
        src/query.cpp
        src/tstatsegmentation.cpp
        src/rawreferences.cpp
        src/rawsegmentation.cpp
        src/rawsignals.cpp
        src/pybindings.cpp)

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
)
add_executable(tests
        tests/testmain.cpp
        tests/tests.h
)

#target_compile_options(collinear PRIVATE -fpermissive -mavx)
target_compile_options(Collinearity PRIVATE -fpermissive)
target_link_libraries(Collinearity slow5 z)

target_compile_options(tests PRIVATE -fpermissive)
target_link_libraries(tests z)

target_link_libraries(_core PUBLIC slow5 z)
install(TARGETS _core DESTINATION pycollinearity)
