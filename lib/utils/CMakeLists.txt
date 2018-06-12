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
  number.cpp
  term.cpp
  timer.cpp
  utils.cpp)
set(nacs_utils_LINKS m pthread)
if(IS_UNIX)
  # For shm_*
  set(nacs_utils_LINKS ${nacs_utils_LINKS} rt)
endif()
set(nacs_utils_llvm_HDRS
  llvm/codegen.h
  llvm/compile.h
  llvm/utils.h)
if(ENABLE_LLVM)
  set(nacs_utils_SRCS ${nacs_utils_SRCS}
    llvm/codegen.cpp
    llvm/compile.cpp
    llvm/execute.cpp
    llvm/memmgr.cpp
    llvm/mergephi-pass.cpp)
endif()

include_directories("${CMAKE_CURRENT_BINARY_DIR}")
if(ENABLE_LLVM)
  include_directories(${LLVM_INCLUDE_DIRS})
  add_definitions(${LLVM_DEFINITIONS})
endif()

add_definitions("-\"DNACS_EXPORT_LIB_utils()=\"")

add_library(nacs-utils SHARED
  ${nacs_utils_SRCS})

if(ENABLE_LLVM)
  # AFAICT llvm_config is pretty broken so just do it ourselves...
  execute_process(COMMAND llvm-config --libs --link-shared OUTPUT_VARIABLE llvm_libs)
  string(STRIP "${llvm_libs}" llvm_libs)
  target_link_libraries(nacs-utils ${llvm_libs})
  if(IS_UNIX)
    target_link_libraries(nacs-utils dl)
  endif()
endif()
target_link_libraries(nacs-utils ${nacs_utils_LINKS})

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
  if(ENABLE_LLVM)
    install(FILES ${nacs_utils_llvm_HDRS}
      DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/nacs-utils/llvm")
  endif()
endif()

install(TARGETS nacs-utils
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}")
