#

set(nacs_bcc_SOURCES nacs_bcc.cpp)
add_executable(nacs-bcc ${nacs_bcc_SOURCES})
target_link_libraries(nacs-bcc nacs-seq)

set(nacs_bcdump_SOURCES nacs_bcdump.cpp)
add_executable(nacs-bcdump ${nacs_bcdump_SOURCES})
target_link_libraries(nacs-bcdump nacs-seq)

set(nacs_cmdlist_SOURCES nacs_cmdlist.cpp)
add_executable(nacs-cmdlist ${nacs_cmdlist_SOURCES})
target_link_libraries(nacs-cmdlist nacs-seq)

set(nacs_seqdump_SOURCES nacs_seqdump.cpp)
add_executable(nacs-seqdump ${nacs_seqdump_SOURCES})
target_link_libraries(nacs-seqdump nacs-seq)

install(TARGETS nacs-bcdump nacs-bcc nacs-cmdlist nacs-seqdump
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}")
