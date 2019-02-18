#

set(test_schedule_SOURCES test_schedule.cpp)
add_executable(test-schedule ${test_schedule_SOURCES})
target_link_libraries(test-schedule nacs-seq)

set(test_cmdlist_parse_SOURCES test_cmdlist_parse.cpp)
add_executable(test-cmdlist_parse ${test_cmdlist_parse_SOURCES})
target_link_libraries(test-cmdlist_parse nacs-seq)
