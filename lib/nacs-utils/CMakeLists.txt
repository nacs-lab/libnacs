#

if(IS_UNIX)
  add_custom_command(OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/nacs_lib_names.h"
    COMMAND "${CMAKE_COMMAND}"
    "-DOUTPUT_FILE=${CMAKE_CURRENT_BINARY_DIR}/nacs_lib_names.h"
    "-DCMAKE_OBJDUMP=${CMAKE_OBJDUMP}"
    "-DNAMES=OPENLIBM"
    "-DOPENLIBM=${OPENLIBM_LIBRARIES}"
    -P "${CMAKE_CURRENT_SOURCE_DIR}/gen_libnames.cmake"
    DEPENDS gen_libnames.cmake
    VERBATIM)
else()
  configure_file(nacs_lib_names.win.h nacs_lib_names.h COPYONLY)
endif()
add_custom_target(nacs-lib-names ALL
  DEPENDS "${CMAKE_CURRENT_BINARY_DIR}/nacs_lib_names.h")
include_directories("${CMAKE_CURRENT_BINARY_DIR}")

set(nacs_utils_HDRS
  container.h
  dlload.h
  errors.h
  fd_utils.h
  features_aarch32.h
  features_aarch64.h
  features_x86.h
  ir.h
  log.h
  macros.h
  mem.h
  number.h
  processor.h
  signals.h
  streams.h
  term.h
  thread.h
  timer.h
  tweezer.h
  utils.h
  wavemeter.h
  zmq_utils.h)
set(nacs_utils_SRCS
  dlload.cpp
  errors.cpp
  fd_utils.cpp
  interp.cpp
  ir.cpp
  ir-interp-shim.S
  log.cpp
  mem.cpp
  number.cpp
  number_simd.S
  processor.cpp
  streams.cpp
  term.cpp
  thread.cpp
  timer.cpp
  tweezer.cpp
  utils.cpp
  wavemeter.cpp
  zmq_utils.cpp)
set(nacs_utils_LINKS m pthread PkgConfig::ZMQ)
set_source_files_properties(number.cpp
  PROPERTIES COMPILE_FLAGS "-ffp-contract=fast")
if(IS_UNIX AND NOT APPLE)
  # For shm_*
  set(nacs_utils_LINKS ${nacs_utils_LINKS} rt)
endif()
# Include the dependencies in the link flags as well since otherwise
# some dependent libraries couldn't be found.
# Ref https://sysprogs.com/w/fixing-rpath-link-issues-with-cross-compilers/
set(nacs_utils_LINKS_FLAGS)
foreach(link IN LISTS nacs_utils_LINKS)
  # Ignore imported target for now since the user of the pc file won't have it
  if(TARGET ${link})
    continue()
  endif()
  set(nacs_utils_LINKS_FLAGS "${nacs_utils_LINKS_FLAGS} -l${link}")
endforeach()
set(nacs_utils_llvm_HDRS
  llvm/analysis.h
  llvm/codegen.h
  llvm/compile.h
  llvm/execute.h
  llvm/mergephi.h
  llvm/passes.h
  llvm/utils.h
  llvm/vectorize.h)
set(nacs_utils_llvm_SRCS
  llvm/analysis.cpp
  llvm/codegen.cpp
  llvm/compile.cpp
  llvm/elim_macho_prefix.cpp
  llvm/execute.cpp
  llvm/global_rename.cpp
  llvm/inst_simplify.cpp
  llvm/lower_vector.cpp
  llvm/memmgr.cpp
  llvm/mergephi.cpp
  llvm/passes.cpp
  llvm/utils.cpp
  llvm/vector_abi.cpp
  llvm/vectorize.cpp)

configure_file(config.h.in config.h @ONLY)
set(nacs_utils_HDRS ${nacs_utils_HDRS} "${CMAKE_CURRENT_BINARY_DIR}/config.h")

if(ENABLE_LLVM)
  set(nacs_utils_SRCS ${nacs_utils_SRCS} ${nacs_utils_llvm_SRCS})
  set_property(SOURCE ${nacs_utils_llvm_SRCS} APPEND_STRING
    PROPERTY COMPILE_FLAGS " ${LLVM_DEFINITIONS} -Wno-conversion")
  # The INCLUDE_DIRECTORIES source property conflicts with the include directories
  # from imported target. Don't use it until it's fixed.
  # See https://gitlab.kitware.com/cmake/cmake/-/issues/22033
  if(TRUE) # (${CMAKE_VERSION} VERSION_LESS "3.11.0")
    include_directories(${LLVM_INCLUDE_DIRS})
  else()
    set_property(SOURCE ${nacs_utils_llvm_SRCS} APPEND
      PROPERTY INCLUDE_DIRECTORIES ${LLVM_INCLUDE_DIRS})
  endif()
else()
  set(nacs_utils_SRCS ${nacs_utils_SRCS} llvm/cpu_nollvm.cpp)
endif()

add_definitions("-\"DNACS_EXPORT_LIB_utils()=\"")

add_library(nacs-utils SHARED
  ${nacs_utils_SRCS})
add_dependencies(nacs-utils nacs-lib-names)
if(BUILD_STATIC_LIB)
  add_library(nacs-utils-static STATIC
    ${nacs_utils_SRCS})
  set_property(TARGET nacs-utils-static APPEND_STRING
    PROPERTY COMPILE_FLAGS " -DNACS_STATIC_ONLY=1")
  add_dependencies(nacs-utils-static nacs-lib-names)
endif()

function(link_nacs_utils)
  target_link_libraries(nacs-utils ${ARGV})
  if(BUILD_STATIC_LIB)
    target_link_libraries(nacs-utils-static ${ARGV})
  endif()
endfunction()

if(ENABLE_LLVM)
  if(STATIC_LLVM)
    # The list below is obtained from trail-and-error.
    # I somehow couldn't find a very clear document on which function
    # belongs to which library...
    set(llvm_libs LLVMCodeGen LLVMExecutionEngine LLVMPasses)
    foreach(target ${LLVM_TARGETS_TO_BUILD})
      if(TARGET "LLVM${target}AsmParser")
        set(llvm_libs ${llvm_libs} "LLVM${target}AsmParser")
      endif()
      set(llvm_libs ${llvm_libs} "LLVM${target}CodeGen")
    endforeach()
  else()
    set(llvm_libs LLVM)
  endif()
  # Link LLVM privately, this makes sure that we are not accidentally linking to LLVM
  # in any other libraries. Additional linking is dangerous for statically linked LLVM
  # since the two libraries/programs will be linking to different copies of LLVM
  # that cannot interact with each other safely.
  # See also comments in `lib/nacs-utils/llvm/utils.h`
  target_link_libraries(nacs-utils PRIVATE ${llvm_libs})
  if(BUILD_STATIC_LIB)
    target_link_libraries(nacs-utils-static PUBLIC ${llvm_libs})
  endif()
  if(MINGW AND STATIC_LLVM)
    link_nacs_utils(PRIVATE ole32 uuid z)
  endif()
  if(ENABLE_SIMD)
    link_nacs_utils(PRIVATE ${SLEEF_LIBRARIES})
  endif()
endif()
if(IS_UNIX)
  link_nacs_utils(PRIVATE dl)
endif()
link_nacs_utils(PUBLIC ${nacs_utils_LINKS})

set_target_properties(nacs-utils PROPERTIES
  VERSION "${MAJOR_VERSION}.${MINOR_VERSION}"
  SOVERSION "${MAJOR_VERSION}"
  COMPILE_FLAGS "-fvisibility=hidden"
  LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")

if(BUILD_STATIC_LIB)
  set_target_properties(nacs-utils-static PROPERTIES
    COMPILE_FLAGS "-fvisibility=hidden"
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")
endif()

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
