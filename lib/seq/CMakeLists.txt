#

set(nacs_seq_HDRS)
set(nacs_seq_zynq_HDRS
  zynq/bytecode.h
  zynq/cmdlist.h
  zynq/seq.h)
set(nacs_seq_SRCS
  zynq/bytecode.cpp
  zynq/cmdlist.cpp
  zynq/parser.cpp
  zynq/seq.cpp)
set(nacs_seq_LINKS)

if(SEQ_IN_UTILS)
  add_library(nacs-seq SHARED anchor.cpp)
  target_link_libraries(nacs-seq PUBLIC nacs-utils)

  foreach(SRC ${nacs_seq_SRCS})
    get_filename_component(SRC "${SRC}" ABSOLUTE)
    target_sources(nacs-utils PRIVATE "${SRC}")
  endforeach()

  target_link_libraries(nacs-utils PUBLIC ${nacs_seq_LINKS})
else()
  add_library(nacs-seq SHARED
    ${nacs_seq_SRCS})

  target_link_libraries(nacs-seq PUBLIC nacs-utils ${nacs_seq_LINKS})
endif()

set_target_properties(nacs-seq PROPERTIES
  VERSION "${MAJOR_VERSION}.${MINOR_VERSION}"
  SOVERSION "${MAJOR_VERSION}"
  COMPILE_FLAGS "-fvisibility=hidden"
  LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")

set_property(SOURCE ${nacs_seq_SRCS} APPEND_STRING
  PROPERTY COMPILE_FLAGS " -\"DNACS_EXPORT_LIB_seq()=\"")

configure_file(nacs-seq.pc.in nacs-seq.pc @ONLY)
if(ENABLE_INSTALL_NONLIB)
  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/nacs-seq.pc
    DESTINATION ${INSTALL_PKGCONFIG_DIR})
  install(FILES ${nacs_seq_HDRS}
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/nacs-seq")
  install(FILES ${nacs_seq_zynq_HDRS}
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/nacs-seq/zynq")
endif()

install(TARGETS nacs-seq
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}")
