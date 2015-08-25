#

set(test_getver_SOURCES test_getver.cpp)
add_executable(test-getver ${test_getver_SOURCES})
target_link_libraries(test-getver nacs-kernel)

set(test_pulse_ctl_SOURCES test_pulse_ctl.cpp)
add_executable(test-pulse_ctl ${test_pulse_ctl_SOURCES})
target_link_libraries(test-pulse_ctl nacs-kernel)
