#

set(nacs_fpga_bcc_SOURCES nacs_fpga_bcc.cpp)
add_executable(nacs-fpga-bcc ${nacs_fpga_bcc_SOURCES})
target_link_libraries(nacs-fpga-bcc nacs-seq)

set(nacs_fpga_bcdump_SOURCES nacs_fpga_bcdump.cpp)
add_executable(nacs-fpga-bcdump ${nacs_fpga_bcdump_SOURCES})
target_link_libraries(nacs-fpga-bcdump nacs-seq)

set(nacs_fpga_cmdlist_SOURCES nacs_fpga_cmdlist.cpp)
add_executable(nacs-fpga-cmdlist ${nacs_fpga_cmdlist_SOURCES})
target_link_libraries(nacs-fpga-cmdlist nacs-seq)

set(nacs_fpga_seqdump_SOURCES nacs_fpga_seqdump.cpp)
add_executable(nacs-fpga-seqdump ${nacs_fpga_seqdump_SOURCES})
target_link_libraries(nacs-fpga-seqdump nacs-seq)

set(nacs_cpudump_SOURCES nacs_cpudump.cpp)
add_executable(nacs-cpudump ${nacs_cpudump_SOURCES})
target_link_libraries(nacs-cpudump nacs-utils)

install(TARGETS nacs-fpga-bcdump nacs-fpga-bcc nacs-fpga-cmdlist nacs-fpga-seqdump nacs-cpudump
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}")
