#

set(test_anyptr_SOURCES test_anyptr.cpp)
add_executable(test-anyptr ${test_anyptr_SOURCES})
target_link_libraries(test-anyptr nacs-utils)

set(test_bits_SOURCES test_bits.cpp)
add_executable(test-bits ${test_bits_SOURCES})
target_link_libraries(test-bits nacs-utils)

set(test_buff_streams_SOURCES test_buff_streams.cpp)
add_executable(test-buff_streams ${test_buff_streams_SOURCES})
target_link_libraries(test-buff_streams nacs-utils)

set(test_const_stream_SOURCES test_const_stream.cpp)
add_executable(test-const_stream ${test_const_stream_SOURCES})
target_link_libraries(test-const_stream nacs-utils)

set(test_filterqueue_SOURCES test_filterqueue.cpp)
add_executable(test-filterqueue ${test_filterqueue_SOURCES})
target_link_libraries(test-filterqueue nacs-utils)

set(test_fixedqueue_SOURCES test_fixedqueue.cpp)
add_executable(test-fixedqueue ${test_fixedqueue_SOURCES})
target_link_libraries(test-fixedqueue nacs-utils)

set(test_interp_SOURCES test_interp.cpp)
add_executable(test-interp ${test_interp_SOURCES})
target_link_libraries(test-interp nacs-utils)

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

set(test_syntax_error_SOURCES test_syntax_error.cpp)
add_executable(test-syntax_error ${test_syntax_error_SOURCES})
target_link_libraries(test-syntax_error nacs-utils)

set(test_term_SOURCES test_term.cpp)
add_executable(test-term ${test_term_SOURCES})
target_link_libraries(test-term nacs-utils)

set(test_thread_SOURCES test_thread.cpp)
add_executable(test-thread ${test_thread_SOURCES})
target_link_libraries(test-thread nacs-utils)

set(test_wavemeter_SOURCES test_wavemeter.cpp)
add_executable(test-wavemeter ${test_wavemeter_SOURCES})
target_link_libraries(test-wavemeter nacs-utils)

set(test_wavemeter2_SOURCES test_wavemeter2.cpp)
add_executable(test-wavemeter2 ${test_wavemeter2_SOURCES})
target_link_libraries(test-wavemeter2 nacs-utils)

if(ENABLE_LLVM)
  set(test_codegen_SOURCES test_codegen.cpp)
  add_executable(test-codegen ${test_codegen_SOURCES})
  target_link_libraries(test-codegen nacs-utils)

  set_property(TARGET test-codegen PROPERTY ${LLVM_INCLUDE_DIRS})
  set_property(SOURCE ${test_codegen_SOURCES} APPEND_STRING
    PROPERTY COMPILE_FLAGS " ${LLVM_DEFINITIONS} -Wno-conversion")
endif()
