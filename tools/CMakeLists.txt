#

set(nacs_zynq_bc_SOURCES nacs_zynq_bc.cpp)
add_executable(nacs-zynq-bc ${nacs_zynq_bc_SOURCES})
target_link_libraries(nacs-zynq-bc nacs-seq)

set(nacs_zynq_cmdlist_SOURCES nacs_zynq_cmdlist.cpp)
add_executable(nacs-zynq-cmdlist ${nacs_zynq_cmdlist_SOURCES})
target_link_libraries(nacs-zynq-cmdlist nacs-seq)

set(nacs_cpudump_SOURCES nacs_cpudump.cpp)
add_executable(nacs-cpudump ${nacs_cpudump_SOURCES})
target_link_libraries(nacs-cpudump nacs-utils)

install(TARGETS nacs-zynq-bc nacs-zynq-cmdlist nacs-cpudump
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}")

if(ENABLE_LLVM)
  set(nacs_seq_tool_SOURCES nacs_seq.cpp)
  add_executable(nacs-seq-tool ${nacs_seq_tool_SOURCES})
  target_link_libraries(nacs-seq-tool nacs-seq)
  set_target_properties(nacs-seq-tool PROPERTIES OUTPUT_NAME nacs-seq)
  set_property(TARGET nacs-seq-tool APPEND
    PROPERTY INCLUDE_DIRECTORIES ${LLVM_INCLUDE_DIRS})
  set_property(TARGET nacs-seq-tool APPEND_STRING
    PROPERTY COMPILE_FLAGS " ${LLVM_DEFINITIONS} -Wno-conversion")
  install(TARGETS nacs-seq-tool
    RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}")
endif()
