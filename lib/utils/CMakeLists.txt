set(nacs_utils_HDRS
  container.h
  errors.h
  fd_utils.h
  ir.h
  log.h
  macros.h
  mem.h
  number.h
  streams.h
  term.h
  thread.h
  timer.h
  utils.h
  wavemeter.h
  zmq_utils.h)
set(nacs_utils_SRCS
  errors.cpp
  fd_utils.cpp
  interp.cpp
  ir.cpp
  ir-interp-shim.S
  log.cpp
  mem.cpp
  number.cpp
  streams.cpp
  term.cpp
  timer.cpp
  utils.cpp
  wavemeter.cpp)
set(nacs_utils_LINKS m pthread)
if(IS_UNIX)
  # For shm_*
  set(nacs_utils_LINKS ${nacs_utils_LINKS} rt)
endif()
# Include the dependencies in the link flags as well since otherwise
# some dependent libraries couldn't be found.
# Ref https://sysprogs.com/w/fixing-rpath-link-issues-with-cross-compilers/
set(nacs_utils_LINKS_FLAGS)
foreach(link IN LISTS nacs_utils_LINKS)
  set(nacs_utils_LINKS_FLAGS "${nacs_utils_LINKS_FLAGS} -l${link}")
endforeach()
set(nacs_utils_llvm_HDRS
  llvm/codegen.h
  llvm/compile.h
  llvm/utils.h)
set(nacs_utils_llvm_SRCS
  llvm/codegen.cpp
  llvm/compile.cpp
  llvm/execute.cpp
  llvm/memmgr.cpp
  llvm/mergephi-pass.cpp)

include_directories("${CMAKE_CURRENT_BINARY_DIR}")
if(ENABLE_LLVM)
  set(nacs_utils_SRCS ${nacs_utils_SRCS} ${nacs_utils_llvm_SRCS})
  include_directories(${LLVM_INCLUDE_DIRS})
  set_property(SOURCE ${nacs_utils_llvm_SRCS} APPEND_STRING
    PROPERTY COMPILE_FLAGS " ${LLVM_DEFINITIONS} -Wno-conversion")
endif()

add_definitions("-\"DNACS_EXPORT_LIB_utils()=\"")

add_library(nacs-utils SHARED
  ${nacs_utils_SRCS})

if(ENABLE_LLVM)
  # AFAICT llvm_config is pretty broken so just do it ourselves...
  if(STATIC_LLVM)
    set(LLVM_CONFIG_FLAGS --link-static)
  else()
    set(LLVM_CONFIG_FLAGS --link-shared)
  endif()
  execute_process(COMMAND llvm-config --libs ${LLVM_CONFIG_FLAGS} OUTPUT_VARIABLE llvm_libs)
  string(STRIP "${llvm_libs}" llvm_libs)
  target_link_libraries(nacs-utils ${llvm_libs})
  if(IS_UNIX)
    target_link_libraries(nacs-utils dl)
  endif()
  if(MINGW AND STATIC_LLVM)
    target_link_libraries(nacs-utils ole32 uuid z)
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
