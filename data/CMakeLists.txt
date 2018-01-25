#

if(ENABLE_WINDOWS_PYWRAPPER)
  configure_file(__init__.py.in __init__.py @ONLY)
  install(FILES "${CMAKE_CURRENT_BINARY_DIR}/__init__.py" DESTINATION "${CMAKE_INSTALL_BINDIR}")
endif()
