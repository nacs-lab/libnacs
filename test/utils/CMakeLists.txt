#

function(add_utils_test target)
  if(${ARGC} EQUAL 1)
    set(SOURCES "test_${target}.cpp")
  else()
    set(SOURCES ${ARGN})
  endif()
  add_executable(test-${target} ${SOURCES})
  target_link_libraries(test-${target} nacs-utils)
endfunction()

add_utils_test(anyptr)
add_utils_test(bits)
add_utils_test(buff_streams)
add_utils_test(const_stream)
add_utils_test(cpu)
add_utils_test(datapipe)
add_utils_test(dl)
add_utils_test(filterqueue)
add_utils_test(fixedqueue)
add_utils_test(interp)

foreach(TESTSET RANGE 8)
  add_utils_test(ir_cc-${TESTSET} test_ir_cc.cpp)
  set_target_properties("test-ir_cc-${TESTSET}" PROPERTIES
    COMPILE_DEFINITIONS "IR_CC_TESTSET=${TESTSET}")
endforeach()

if(CMAKE_SYSTEM_NAME STREQUAL Linux)
  add_utils_test(phy_addr)
endif()

add_utils_test(small_alloc)
add_utils_test(syntax_error)
add_utils_test(term)
add_utils_test(thread)
add_utils_test(wavemeter)
add_utils_test(wavemeter2)

if(ENABLE_LLVM)
  function(add_llvm_test target)
    add_utils_test(${ARGV})
    set_property(TARGET test-${target} APPEND
      PROPERTY INCLUDE_DIRECTORIES ${LLVM_INCLUDE_DIRS})
    set_property(TARGET test-${target} APPEND_STRING
      PROPERTY COMPILE_FLAGS " ${LLVM_DEFINITIONS} -Wno-conversion")
  endfunction()

  add_llvm_test(codegen)
  add_llvm_test(codegen_dump)
  add_llvm_test(vec_codegen)
endif()
