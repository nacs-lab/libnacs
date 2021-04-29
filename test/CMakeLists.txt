add_definitions(-UNDEBUG)

find_package(Catch2 CONFIG)
if(NOT TARGET Catch2::Catch2)
  add_library(Catch2::Catch2 INTERFACE IMPORTED)
  set_target_properties(Catch2::Catch2 PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${CMAKE_CURRENT_SOURCE_DIR}/bundle")
endif()

configure_file(winpath_helper.h.in winpath_helper.h @ONLY)
include_directories("${CMAKE_CURRENT_BINARY_DIR}")

add_subdirectory(utils)
add_subdirectory(seq)
if(ENABLE_KNACS)
  add_subdirectory(kernel)
endif()
