cmake_minimum_required(VERSION 3.0.0)
project(elliptic_curves VERSION 0.1.0)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_EXPORT_COMPILE_COMMANDS true)
set(CMAKE_CXX_FLAGS "-Ofast -flto -ffast-math -Wall -Wextra -Wno-unused-parameter -pedantic")
set(CMAKE_EXE_LINKER_FLAGS "-Ofast -flto -ffast-math")

include_directories(elliptic_curves PUBLIC . base)

add_executable(elliptic_curves main.cpp)
# set_target_properties(elliptic_curves PROPERTIES COMPILE_FLAGS "-DNDEBUG")
