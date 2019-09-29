#

set(nacs_zynq_bcc_SOURCES nacs_zynq_bcc.cpp)
add_executable(nacs-zynq-bcc ${nacs_zynq_bcc_SOURCES})
target_link_libraries(nacs-zynq-bcc nacs-seq)

set(nacs_zynq_bcdump_SOURCES nacs_zynq_bcdump.cpp)
add_executable(nacs-zynq-bcdump ${nacs_zynq_bcdump_SOURCES})
target_link_libraries(nacs-zynq-bcdump nacs-seq)

set(nacs_zynq_cmdlist_SOURCES nacs_zynq_cmdlist.cpp)
add_executable(nacs-zynq-cmdlist ${nacs_zynq_cmdlist_SOURCES})
target_link_libraries(nacs-zynq-cmdlist nacs-seq)

set(nacs_zynq_seqdump_SOURCES nacs_zynq_seqdump.cpp)
add_executable(nacs-zynq-seqdump ${nacs_zynq_seqdump_SOURCES})
target_link_libraries(nacs-zynq-seqdump nacs-seq)

set(nacs_cpudump_SOURCES nacs_cpudump.cpp)
add_executable(nacs-cpudump ${nacs_cpudump_SOURCES})
target_link_libraries(nacs-cpudump nacs-utils)

install(TARGETS nacs-zynq-bcdump nacs-zynq-bcc nacs-zynq-cmdlist nacs-zynq-seqdump nacs-cpudump
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}")
