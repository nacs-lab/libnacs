#

set(test_schedule_SOURCES test_schedule.cpp)
add_executable(test-schedule ${test_schedule_SOURCES})
target_link_libraries(test-schedule nacs-seq)
