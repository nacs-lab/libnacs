#

include(PythonMacros)

if(IS_WINDOWS)
  cmake_utils_cygpath_m(INSTALL_WINDOWS_BINDIR "${CMAKE_INSTALL_FULL_BINDIR}")
  configure_file(handle_win.py.in handle.py @ONLY)
else()
  configure_file(handle_lin.py.in handle.py @ONLY)
endif()

python_install_module(libnacs "${CMAKE_CURRENT_BINARY_DIR}/handle.py"
  __init__.py
  wavemeter.py)
