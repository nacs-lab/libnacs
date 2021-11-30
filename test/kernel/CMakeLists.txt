#

function(add_kernel_test target)
  if(${ARGC} EQUAL 1)
    set(SOURCES "test_${target}.cpp")
  else()
    set(SOURCES ${ARGN})
  endif()
  add_executable(test-${target} ${SOURCES})
  target_link_libraries(test-${target} nacs-kernel Catch2::Catch2)
endfunction()

add_kernel_test(getver)
add_kernel_test(pulse_ctl)
add_kernel_test(dma_buff)
