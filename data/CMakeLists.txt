#

if(ENABLE_WINDOWS_PYWRAPPER)
  execute_process(COMMAND cygpath -m "${CMAKE_INSTALL_FULL_BINDIR}"
    OUTPUT_VARIABLE INSTALL_WINDOWS_BINDIR)
  string(STRIP "${INSTALL_WINDOWS_BINDIR}" INSTALL_WINDOWS_BINDIR)
  configure_file(__init__.py.in __init__.py @ONLY)
  install(FILES "${CMAKE_CURRENT_BINARY_DIR}/__init__.py"
    DESTINATION "${INSTALL_PYTHON_SITE}/libnacs")
endif()
