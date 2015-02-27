set(test_thread_SOURCES test_thread.cpp)

add_executable(test-thread ${test_thread_SOURCES})

target_link_libraries(test-thread nacs-utils)
