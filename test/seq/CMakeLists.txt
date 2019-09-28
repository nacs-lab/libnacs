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
