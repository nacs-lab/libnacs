set(nacs_utils_HDRS
  container.h
  fd_utils.h
  ir.h
  log.h
  macros.h
  mem.h
  number.h
  term.h
  thread.h
  timer.h
  utils.h
  zmq_utils.h)
set(nacs_utils_SRCS
  fd_utils.cpp
  interp.cpp
  ir.cpp
  ir-interp-shim.S
  log.cpp
  mem.cpp
  term.cpp
  timer.cpp)
set(nacs_utils_LINKS m pthread)

include_directories("${CMAKE_CURRENT_BINARY_DIR}")
if(ENABLE_LLVM)
  include_directories(${LLVM_INCLUDE_DIRS})
  add_definitions(${LLVM_DEFINITIONS})
endif()

add_definitions("-\"DNACS_EXPORT_LIB_utils()=\"")

add_library(nacs-utils SHARED
  ${nacs_utils_SRCS})

if(ENABLE_LLVM)
  llvm_config(nacs-utils USE_SHARED all)
endif()
target_link_libraries(nacs-utils PUBLIC ${nacs_utils_LINKS})

set_target_properties(nacs-utils PROPERTIES
  VERSION "${MAJOR_VERSION}.${MINOR_VERSION}"
  SOVERSION "${MAJOR_VERSION}"
  COMPILE_FLAGS "-fvisibility=hidden"
  LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")

configure_file(nacs-utils.pc.in nacs-utils.pc @ONLY)
if(ENABLE_INSTALL_NONLIB)
  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/nacs-utils.pc
    DESTINATION ${INSTALL_PKGCONFIG_DIR})
  install(FILES ${nacs_utils_HDRS}
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/nacs-utils")
endif()

install(TARGETS nacs-utils
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}")
