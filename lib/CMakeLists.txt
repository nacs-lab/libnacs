set(nacs_utils_HDRS
  utils.h
  macros.h)

install(FILES ${nacs_utils_HDRS}
  DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/nacs-utils")
