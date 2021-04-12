#

function(add_seq_test target)
  cmake_parse_arguments(SEQ_TEST "NOTEST;EXCLUDE" "" "" ${ARGN})
  if(SEQ_TEST_UNPARSED_ARGUMENTS)
    set(SOURCES ${SEQ_TEST_UNPARSED_ARGUMENTS})
  else()
    set(SOURCES "test_${target}.cpp")
  endif()
  add_executable(test-${target} ${SOURCES})
  target_link_libraries(test-${target} nacs-seq Catch2::Catch2)
  if(NOT SEQ_TEST_NOTEST)
    if(SEQ_TEST_EXCLUDE)
      add_test(NAME test/excluded/seq/${target} COMMAND test-${target})
    else()
      add_test(NAME test/normal/seq/${target} COMMAND test-${target})
    endif()
  endif()
endfunction()

add_seq_test(schedule NOTEST)
add_seq_test(cmdlist_parse NOTEST)
add_test(NAME test/normal/seq/cmdlist_parse
  COMMAND ${CMAKE_COMMAND} -E env "TEST_SOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR}"
  "$<TARGET_FILE:test-cmdlist_parse>")
add_seq_test(host_seq)
add_seq_test(host_seq_error)
add_seq_test(nidaq_data_gen)
add_seq_test(zynq_bc_gen)

if(ENABLE_LLVM)
  function(add_llvm_test target)
    add_seq_test(${ARGV})
    set_property(TARGET test-${target} APPEND
      PROPERTY INCLUDE_DIRECTORIES ${LLVM_INCLUDE_DIRS})
    set_property(TARGET test-${target} APPEND_STRING
      PROPERTY COMPILE_FLAGS " ${LLVM_DEFINITIONS} -Wno-conversion")
  endfunction()

  add_llvm_test(env)
  add_llvm_test(optimize)

  add_llvm_test(event_time)
  add_llvm_test(pulse)
  add_llvm_test(seq)
  add_llvm_test(seq_check)
  add_llvm_test(seq_optimize)

  add_llvm_test(builder)
  add_llvm_test(builder_check)
  add_llvm_test(builder_deserial)

  add_llvm_test(compiler)
  add_llvm_test(zynq_backend)
endif()
