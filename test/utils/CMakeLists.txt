set(test_thread_SOURCES test_thread.cpp)
add_executable(test-thread ${test_thread_SOURCES})
target_link_libraries(test-thread nacs-utils)

set(test_bits_SOURCES test_bits.cpp)
add_executable(test-bits ${test_bits_SOURCES})
target_link_libraries(test-bits nacs-utils)

set(test_fifo_SOURCES test_fifo.cpp)
add_executable(test-fifo ${test_fifo_SOURCES})
target_link_libraries(test-fifo nacs-utils)

set(test_phy_addr_SOURCES test_phy_addr.cpp)
add_executable(test-phy_addr ${test_phy_addr_SOURCES})
target_link_libraries(test-phy_addr nacs-utils)
