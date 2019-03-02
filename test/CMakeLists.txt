add_definitions(-UNDEBUG)
configure_file(winpath_helper.h.in winpath_helper.h @ONLY)
include_directories("${CMAKE_CURRENT_BINARY_DIR}")

add_subdirectory(utils)
add_subdirectory(seq)
if(ENABLE_KNACS)
  add_subdirectory(kernel)
endif()
