#

set(nacs_bcdump_SOURCES nacs_bcdump.cpp)
add_executable(nacs-bcdump ${nacs_bcdump_SOURCES})
target_link_libraries(nacs-bcdump nacs-seq)

install(TARGETS nacs-bcdump
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}")
