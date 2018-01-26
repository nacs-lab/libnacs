#

include(PythonMacros)

if(IS_WINDOWS)
  execute_process(COMMAND cygpath -m "${CMAKE_INSTALL_FULL_BINDIR}"
    OUTPUT_VARIABLE INSTALL_WINDOWS_BINDIR)
  string(STRIP "${INSTALL_WINDOWS_BINDIR}" INSTALL_WINDOWS_BINDIR)
  configure_file(handle_win.py.in handle.py @ONLY)
else()
  configure_file(handle_lin.py.in handle.py @ONLY)
endif()

python_install_module(libnacs "${CMAKE_CURRENT_BINARY_DIR}/handle.py" __init__.py)
