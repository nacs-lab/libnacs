#

include(PythonMacros)

if(IS_WINDOWS)
  configure_file(handle_win.py.in handle.py @ONLY)
else()
  configure_file(handle_lin.py.in handle.py @ONLY)
endif()

python_install_module(libnacs "${CMAKE_CURRENT_BINARY_DIR}/handle.py"
  __init__.py
  tweezer.py
  wavemeter.py)
