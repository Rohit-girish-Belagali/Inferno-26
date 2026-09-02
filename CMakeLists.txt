cmake_minimum_required(VERSION 3.20)
project(inferno_solver LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

add_compile_options(-Wall -Wextra)

include_directories(${CMAKE_SOURCE_DIR})

add_library(inferno_core
  io/mps_reader.cpp
  io/solution_writer.cpp
  checker/checker.cpp
  simplex/dense_simplex.cpp
)
target_include_directories(inferno_core PUBLIC ${CMAKE_SOURCE_DIR})

add_executable(solver src/main.cpp)
target_link_libraries(solver PRIVATE inferno_core)

enable_testing()
add_subdirectory(tests)
