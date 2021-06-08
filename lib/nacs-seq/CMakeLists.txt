#

set(nacs_seq_HDRS
  error.h
  host_seq.h
  sign.h)
set(nacs_seq_zynq_HDRS
  zynq/backend.h
  zynq/bc_gen.h
  zynq/bytecode.h
  zynq/cmdlist.h
  zynq/pulse_time.h
  zynq/legacy_seq.h)
set(nacs_seq_nidaq_HDRS
  nidaq/backend.h
  nidaq/data_gen.h)
set(nacs_seq_awg_HDRS
  awg/backend.h)
set(nacs_seq_SRCS
  error.cpp
  host_seq.cpp
  nidaq/data_gen.cpp
  zynq/bc_gen.cpp
  zynq/bytecode.cpp
  zynq/cmdlist.cpp
  zynq/parser.cpp
  zynq/legacy_seq.cpp)
set(nacs_seq_LINKS)
set(nacs_seq_llvm_HDRS
  basic_seq.h
  builder.h
  compiler.h
  device.h
  env.h
  event_time.h
  manager.h
  pulse.h
  seq.h)
set(nacs_seq_llvm_SRCS
  awg/backend.cpp
  basic_seq.cpp
  builder.cpp
  compiler.cpp
  device.cpp
  env.cpp
  event_time.cpp
  manager.cpp
  nidaq/backend.cpp
  pulse.cpp
  seq.cpp
  var.cpp
  zynq/backend.cpp)
set_source_files_properties(zynq/bc_gen.cpp
  PROPERTIES COMPILE_FLAGS "-ffp-contract=fast")

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
    target_link_libraries(nacs-utils PRIVATE yaml-cpp PkgConfig::ZMQ)
    if(NOT APPLE)
      target_link_libraries(nacs-utils PRIVATE tbb)
    endif()
  endif()
else()
  add_library(nacs-seq SHARED
    ${nacs_seq_SRCS})

  target_link_libraries(nacs-seq PUBLIC nacs-utils ${nacs_seq_LINKS})

  if(ENABLE_LLVM)
    set_property(SOURCE ${nacs_seq_llvm_SRCS} APPEND_STRING
      PROPERTY COMPILE_FLAGS " ${LLVM_DEFINITIONS} -Wno-conversion")
    # The INCLUDE_DIRECTORIES source property conflicts with the include directories
    # from imported target. Don't use it until it's fixed.
    # See https://gitlab.kitware.com/cmake/cmake/-/issues/22033
    if(TRUE) # (${CMAKE_VERSION} VERSION_LESS "3.11.0")
      include_directories(${LLVM_INCLUDE_DIRS})
    else()
      set_property(SOURCE ${nacs_seq_llvm_SRCS} APPEND
        PROPERTY INCLUDE_DIRECTORIES ${LLVM_INCLUDE_DIRS})
    endif()
    target_link_libraries(nacs-seq PRIVATE LLVM yaml-cpp PkgConfig::ZMQ)
    if(NOT APPLE)
      target_link_libraries(nacs-seq PRIVATE tbb)
    endif()
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
  install(FILES ${nacs_seq_nidaq_HDRS}
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/nacs-seq/nidaq")
  install(FILES ${nacs_seq_awg_HDRS}
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/nacs-seq/awg")
endif()

install(TARGETS nacs-seq
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}")
