set(nacs_seq_HDRS
  bytecode.h
  cmdlist.h
  misc.h
  seq.h)
set(nacs_seq_SRCS
  bytecode.cpp
  cmdlist.cpp
  misc.cpp
  parser.cpp
  seq.cpp)
set(nacs_seq_LINKS nacs-utils)

include_directories("${CMAKE_CURRENT_BINARY_DIR}")
add_definitions("-\"DNACS_EXPORT_LIB_seq()=\"")

add_library(nacs-seq SHARED
  ${nacs_seq_SRCS})

target_link_libraries(nacs-seq ${nacs_seq_LINKS})

set_target_properties(nacs-seq PROPERTIES
  VERSION "${MAJOR_VERSION}.${MINOR_VERSION}"
  SOVERSION "${MAJOR_VERSION}"
  COMPILE_FLAGS "-fvisibility=hidden"
  LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")

configure_file(nacs-seq.pc.in nacs-seq.pc @ONLY)
if(ENABLE_INSTALL_NONLIB)
  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/nacs-seq.pc
    DESTINATION ${INSTALL_PKGCONFIG_DIR})
  install(FILES ${nacs_seq_HDRS}
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/nacs-seq")
endif()

install(TARGETS nacs-seq
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}")
