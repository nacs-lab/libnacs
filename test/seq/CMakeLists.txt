#

function(add_seq_test target)
  # COST is a rough estimate of the run time (in seconds on the CI)
  # so that the long running tests are started first.
  cmake_parse_arguments(SEQ_TEST "NOTEST;EXCLUDE;STATIC" "COST" "" ${ARGN})
  if(NOT DEFINED SEQ_TEST_COST)
    set(SEQ_TEST_COST 0.1)
  endif()
  if(SEQ_TEST_UNPARSED_ARGUMENTS)
    set(SOURCES ${SEQ_TEST_UNPARSED_ARGUMENTS})
  else()
    set(SOURCES "test_${target}.cpp")
  endif()
  add_executable(test-${target} ${SOURCES})
  if(SEQ_TEST_STATIC AND BUILD_STATIC_LIB)
    target_link_libraries(test-${target}
      nacs-seq-static nacs-catch-main)
    set_property(TARGET test-${target} APPEND_STRING
      PROPERTY COMPILE_FLAGS " -DNACS_STATIC_ONLY=1")
  else()
    target_link_libraries(test-${target} nacs-seq nacs-catch-main)
  endif()
  if(NOT SEQ_TEST_NOTEST)
    if(SEQ_TEST_EXCLUDE)
      set(test_name test/excluded/seq/${target})
    else()
      set(test_name test/normal/seq/${target})
    endif()
    add_test(NAME ${test_name} COMMAND test-${target})
    set_tests_properties(${test_name} PROPERTIES COST ${SEQ_TEST_COST})
  endif()
endfunction()

add_seq_test(schedule NOTEST)
add_seq_test(cmdlist_parse NOTEST)
add_test(NAME test/normal/seq/cmdlist_parse
  COMMAND ${CMAKE_COMMAND} -E env "TEST_SOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR}"
  "$<TARGET_FILE:test-cmdlist_parse>")
set_tests_properties(test/normal/seq/cmdlist_parse PROPERTIES COST 0.1)
add_seq_test(host_seq)
add_seq_test(host_seq_error)
add_seq_test(nidaq_data_gen)
add_seq_test(zynq_bc_gen COST 1)

if(ENABLE_LLVM)
  function(add_llvm_test target)
    add_seq_test(${ARGV} STATIC)
    set_property(TARGET test-${target} APPEND
      PROPERTY INCLUDE_DIRECTORIES ${LLVM_INCLUDE_DIRS})
    set_property(TARGET test-${target} APPEND_STRING
      PROPERTY COMPILE_FLAGS " ${LLVM_DEFINITIONS} -Wno-conversion")
    if(NOT BUILD_STATIC_LIB)
      target_link_libraries(test-${target} LLVM)
    endif()
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
  add_llvm_test(nidaq_backend)
  add_llvm_test(zynq_backend)
endif()
