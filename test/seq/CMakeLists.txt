#

set(test_schedule_SOURCES test_schedule.cpp)
add_executable(test-schedule ${test_schedule_SOURCES})
target_link_libraries(test-schedule nacs-seq)

set(test_schedule2_SOURCES test_schedule2.cpp)
add_executable(test-schedule2 ${test_schedule2_SOURCES})
target_link_libraries(test-schedule2 nacs-seq)
