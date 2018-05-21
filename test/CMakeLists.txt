
include_directories(
    ${EIGEN3_INCLUDE_DIRS}
    )



set(ProMP_tests
    Phase_test.cpp
    ProMP_test_main.cpp)



add_executable (ProMP_cpp_test ${ProMP_tests})
target_link_libraries(ProMP_cpp_test ProMP)

enable_testing()
add_test(ProMP_test ProMP_cpp_test)
