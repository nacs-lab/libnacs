#

set(test_getver_SOURCES test_getver.cpp)
add_executable(test-getver ${test_getver_SOURCES})
target_link_libraries(test-getver nacs-kernel)
