add_subdirectory(nacs-utils)
add_subdirectory(nacs-seq)
if(ENABLE_KNACS)
  add_subdirectory(nacs-kernel)
endif()
