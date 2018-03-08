#

add_definitions(-UNDEBUG)

set(test_anyptr_SOURCES test_anyptr.cpp)
add_executable(test-anyptr ${test_anyptr_SOURCES})
target_link_libraries(test-anyptr nacs-utils)

set(test_bits_SOURCES test_bits.cpp)
add_executable(test-bits ${test_bits_SOURCES})
target_link_libraries(test-bits nacs-utils)

set(test_fifo_SOURCES test_fifo.cpp)
add_executable(test-fifo ${test_fifo_SOURCES})
target_link_libraries(test-fifo nacs-utils)

set(test_ir_SOURCES test_ir.cpp)
add_executable(test-ir ${test_ir_SOURCES})
target_link_libraries(test-ir nacs-utils)

set(test_ir_cc_SOURCES test_ir_cc.cpp)
foreach(TESTSET RANGE 8)
  add_executable("test-ir_cc-${TESTSET}" ${test_ir_cc_SOURCES})
  target_link_libraries("test-ir_cc-${TESTSET}" nacs-utils)
  set_target_properties("test-ir_cc-${TESTSET}" PROPERTIES
    COMPILE_DEFINITIONS "IR_CC_TESTSET=${TESTSET}")
endforeach()

if(CMAKE_SYSTEM_NAME STREQUAL Linux)
  set(test_phy_addr_SOURCES test_phy_addr.cpp)
  add_executable(test-phy_addr ${test_phy_addr_SOURCES})
  target_link_libraries(test-phy_addr nacs-utils)
endif()

set(test_small_alloc_SOURCES test_small_alloc.cpp)
add_executable(test-small_alloc ${test_small_alloc_SOURCES})
target_link_libraries(test-small_alloc nacs-utils)

set(test_term_SOURCES test_term.cpp)
add_executable(test-term ${test_term_SOURCES})
target_link_libraries(test-term nacs-utils)

set(test_thread_SOURCES test_thread.cpp)
add_executable(test-thread ${test_thread_SOURCES})
target_link_libraries(test-thread nacs-utils)
