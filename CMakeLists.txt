cmake_minimum_required(VERSION 3.0.0)
project(finite_fields VERSION 0.1.0)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_EXPORT_COMPILE_COMMANDS true)
set(CMAKE_CXX_FLAGS "-Ofast -flto -ffast-math -Wall -Wextra -Wno-unused-parameter -pedantic")
set(CMAKE_EXE_LINKER_FLAGS "-Ofast -flto -ffast-math")

set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR})
find_package(GMP REQUIRED)

include_directories(finite_fields PUBLIC ${CMAKE_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/base ${GMP_INCLUDE_DIRS})

add_executable(finite_fields main.cpp)
target_link_libraries(finite_fields gmp gmpxx)
# set_target_properties(finite_fields PROPERTIES COMPILE_FLAGS "-DNDEBUG")
