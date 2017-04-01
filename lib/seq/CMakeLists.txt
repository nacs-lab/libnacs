set(nacs_seq_HDRS
  seq.h
  timing.h
  pulser.h)
set(nacs_seq_SRCS
  timing.cpp
  pulser.cpp)
set(nacs_seq_LINKS nacs-utils)

include_directories("${CMAKE_CURRENT_BINARY_DIR}")

add_library(nacs-seq SHARED
  ${nacs_seq_SRCS})

target_link_libraries(nacs-seq ${nacs_seq_LINKS})

set_target_properties(nacs-seq PROPERTIES
  VERSION 1.1
  SOVERSION 1
  COMPILE_FLAGS "-fvisibility=hidden"
  LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")

install(TARGETS nacs-seq LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}")
install(FILES ${nacs_seq_HDRS}
  DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/nacs-seq")
