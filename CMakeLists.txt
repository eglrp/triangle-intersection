project(triangle_intersection)
cmake_minimum_required(VERSION 2.8)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/build")

FILE(GLOB_RECURSE SRC_LIST "*.cpp")
include_directories("include")

add_library(${PROJECT_NAME} SHARED ${SRC_LIST})

if(UNIX)
    target_compile_options(${PROJECT_NAME} PRIVATE -Wall -Wno-reorder -O3 -m64 -mfpmath=sse -std=c++11 -ftree-vectorizer-verbose=2)
endif(UNIX)
