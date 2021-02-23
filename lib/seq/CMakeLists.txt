#

set(nacs_seq_HDRS)
set(nacs_seq_zynq_HDRS
  zynq/bytecode.h
  zynq/cmdlist.h
  zynq/legacy_seq.h)
set(nacs_seq_SRCS
  zynq/bytecode.cpp
  zynq/cmdlist.cpp
  zynq/parser.cpp
  zynq/legacy_seq.cpp)
set(nacs_seq_LINKS)
set(nacs_seq_llvm_HDRS
  basic_seq.h
  env.h
  error.h
  event_time.h
  pulse.h
  seq.h)
set(nacs_seq_llvm_SRCS
  basic_seq.cpp
  env.cpp
  error.cpp
  event_time.cpp
  pulse.cpp
  seq.cpp
  var.cpp)

if(ENABLE_LLVM)
  set(nacs_seq_SRCS ${nacs_seq_SRCS} ${nacs_seq_llvm_SRCS})
  set(nacs_seq_HDRS ${nacs_seq_HDRS} ${nacs_seq_llvm_HDRS})
endif()

if(SEQ_IN_UTILS)
  add_library(nacs-seq SHARED anchor.cpp)
  target_link_libraries(nacs-seq PUBLIC nacs-utils)

  foreach(SRC ${nacs_seq_SRCS})
    get_filename_component(SRC "${SRC}" ABSOLUTE)
    target_sources(nacs-utils PRIVATE "${SRC}")
  endforeach()

  target_link_libraries(nacs-utils PUBLIC ${nacs_seq_LINKS})

  if(ENABLE_LLVM)
    set_property(TARGET nacs-utils APPEND_STRING
      PROPERTY COMPILE_FLAGS " ${LLVM_DEFINITIONS} -Wno-conversion")
    set_property(TARGET nacs-utils APPEND_STRING
      PROPERTY INCLUDE_DIRECTORIES ${LLVM_INCLUDE_DIRS})
  endif()
else()
  add_library(nacs-seq SHARED
    ${nacs_seq_SRCS})

  target_link_libraries(nacs-seq PUBLIC nacs-utils ${nacs_seq_LINKS})

  if(ENABLE_LLVM)
    set_property(SOURCE ${nacs_seq_llvm_SRCS} APPEND_STRING
      PROPERTY COMPILE_FLAGS " ${LLVM_DEFINITIONS} -Wno-conversion")
    if(${CMAKE_VERSION} VERSION_LESS "3.11.0")
      include_directories(${LLVM_INCLUDE_DIRS})
    else()
      set_property(SOURCE ${nacs_seq_llvm_SRCS} APPEND
        PROPERTY INCLUDE_DIRECTORIES ${LLVM_INCLUDE_DIRS})
    endif()
    target_link_libraries(nacs-seq PRIVATE LLVM)
  endif()
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
