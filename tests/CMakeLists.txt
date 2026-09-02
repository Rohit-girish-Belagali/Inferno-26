add_executable(smoke_test smoke_test.cpp)
target_link_libraries(smoke_test PRIVATE inferno_core)
add_test(NAME smoke_test COMMAND smoke_test)

add_executable(la_test la_test.cpp)
target_link_libraries(la_test PRIVATE inferno_core)
add_test(NAME la_test COMMAND la_test)
