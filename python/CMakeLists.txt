#

include(PythonMacros)

set(NACS_UTILS nacs-utils)
if(SEQ_IN_UTILS)
  set(NACS_SEQ nacs-utils)
else()
  set(NACS_SEQ nacs-seq)
endif()

if(IS_WINDOWS)
  configure_file(handle_win.py.in handle.py @ONLY)
else()
  configure_file(handle_lin.py.in handle.py @ONLY)
endif()

python_install_module(libnacs "${CMAKE_CURRENT_BINARY_DIR}/handle.py"
  __init__.py
  tweezer.py
  msquared.py
  wavemeter.py)
