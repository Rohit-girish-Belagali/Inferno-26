add_executable(smoke_test smoke_test.cpp)
target_link_libraries(smoke_test PRIVATE inferno_core)
add_test(NAME smoke_test COMMAND smoke_test)

add_executable(la_test la_test.cpp)
target_link_libraries(la_test PRIVATE inferno_core)
add_test(NAME la_test COMMAND la_test)

add_executable(revised_simplex_test revised_simplex_test.cpp)
target_link_libraries(revised_simplex_test PRIVATE inferno_core)
add_test(NAME revised_simplex_test COMMAND revised_simplex_test)
