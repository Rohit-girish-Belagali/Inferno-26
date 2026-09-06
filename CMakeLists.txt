cmake_minimum_required(VERSION 3.20)
project(inferno_solver LANGUAGES CXX C)

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
  simplex/revised_simplex.cpp
  simplex/dual_simplex.cpp
  simplex/solve_manager.cpp
  presolve/presolve.cpp
  firstorder/pdlp.cpp
  models/crude_blending.cpp
  models/breadth_models.cpp
  models/industrial_models.cpp
  mip/branch_and_bound.cpp
  qp/admm.cpp
  api/inferno.cpp
  la/scaling.cpp
  la/lu_factors.cpp
  la/markowitz_lu.cpp
  la/lu_solve.cpp
  la/basis_factorization.cpp
)
target_include_directories(inferno_core PUBLIC ${CMAKE_SOURCE_DIR})

# Shared library exposing the C ABI, for the Python bindings and any other
# out-of-process consumer. Built from the same sources as the static
# library rather than wrapping it, since linking a static archive into a
# shared object needs position-independent code throughout and this is
# simpler than arranging that.
add_library(inferno_shared SHARED
  io/mps_reader.cpp
  io/solution_writer.cpp
  checker/checker.cpp
  simplex/dense_simplex.cpp
  simplex/revised_simplex.cpp
  simplex/dual_simplex.cpp
  simplex/solve_manager.cpp
  presolve/presolve.cpp
  firstorder/pdlp.cpp
  models/crude_blending.cpp
  models/breadth_models.cpp
  models/industrial_models.cpp
  mip/branch_and_bound.cpp
  qp/admm.cpp
  api/inferno.cpp
  la/scaling.cpp
  la/lu_factors.cpp
  la/markowitz_lu.cpp
  la/lu_solve.cpp
  la/basis_factorization.cpp
)
target_include_directories(inferno_shared PUBLIC ${CMAKE_SOURCE_DIR})
set_target_properties(inferno_shared PROPERTIES OUTPUT_NAME inferno)

add_executable(solver src/main.cpp)
target_link_libraries(solver PRIVATE inferno_core)

add_executable(refinery_demo src/refinery_demo.cpp)
target_link_libraries(refinery_demo PRIVATE inferno_core)

add_executable(battery src/battery.cpp)
target_link_libraries(battery PRIVATE inferno_core)

enable_testing()
add_subdirectory(tests)
