add_executable(smoke_test smoke_test.cpp)
target_link_libraries(smoke_test PRIVATE inferno_core)
add_test(NAME smoke_test COMMAND smoke_test)

add_executable(la_test la_test.cpp)
target_link_libraries(la_test PRIVATE inferno_core)
add_test(NAME la_test COMMAND la_test)

add_executable(revised_simplex_test revised_simplex_test.cpp)
target_link_libraries(revised_simplex_test PRIVATE inferno_core)
add_test(NAME revised_simplex_test COMMAND revised_simplex_test)

add_executable(presolve_test presolve_test.cpp)
target_link_libraries(presolve_test PRIVATE inferno_core)
add_test(NAME presolve_test COMMAND presolve_test WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})

add_executable(dual_simplex_test dual_simplex_test.cpp)
target_link_libraries(dual_simplex_test PRIVATE inferno_core)
add_test(NAME dual_simplex_test COMMAND dual_simplex_test WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})

add_executable(refinery_test refinery_test.cpp)
target_link_libraries(refinery_test PRIVATE inferno_core)
add_test(NAME refinery_test COMMAND refinery_test WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})

add_executable(qp_test qp_test.cpp)
target_link_libraries(qp_test PRIVATE inferno_core)
add_test(NAME qp_test COMMAND qp_test WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})

# Compiled as C on purpose: the ABI exists so C consumers can use it.
add_executable(api_test api_test.c)
target_link_libraries(api_test PRIVATE inferno_core)
set_target_properties(api_test PROPERTIES LINKER_LANGUAGE CXX)
add_test(NAME api_test COMMAND api_test WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})

# Bindings that are not run are bindings that rot.
add_test(NAME python_bindings_test
         COMMAND python3 ${CMAKE_SOURCE_DIR}/python/test_bindings.py
         WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})

add_executable(breadth_test breadth_test.cpp)
target_link_libraries(breadth_test PRIVATE inferno_core)
add_test(NAME breadth_test COMMAND breadth_test WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
