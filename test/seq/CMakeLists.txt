#

function(add_seq_test target)
  if(${ARGC} EQUAL 1)
    set(SOURCES "test_${target}.cpp")
  else()
    set(SOURCES ${ARGN})
  endif()
  add_executable(test-${target} ${SOURCES})
  target_link_libraries(test-${target} nacs-seq)
endfunction()

add_seq_test(schedule)
add_seq_test(cmdlist_parse)

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
endif()
