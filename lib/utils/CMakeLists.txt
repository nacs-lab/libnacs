set(nacs_utils_HDRS
  utils.h
  macros.h
  log.h
  number.h
  fd_utils.h
  timer.h
  mem.h
  thread.h
  container.h
  ir.h)
set(nacs_utils_SRCS
  log.cpp
  fd_utils.cpp
  timer.cpp
  mem.cpp
  ir.cpp)
set(nacs_utils_LINKS m rt pthread)

include_directories("${CMAKE_CURRENT_BINARY_DIR}")

add_library(nacs-utils SHARED
  ${nacs_utils_SRCS})

target_link_libraries(nacs-utils ${nacs_utils_LINKS})

set_target_properties(nacs-utils PROPERTIES
  VERSION 0.1
  SOVERSION 0
  COMPILE_FLAGS "-fvisibility=hidden"
  LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")

install(TARGETS nacs-utils LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}")
install(FILES ${nacs_utils_HDRS}
  DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/nacs-utils")
