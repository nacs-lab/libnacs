set(nacs_kernel_HDRS
  device.h
  devctl.h
  dma_buffer.h)
set(nacs_kernel_SRCS
  device.cpp
  devctl.cpp
  dma_buffer.cpp)
set(nacs_kernel_LINKS nacs-utils)

include_directories("${CMAKE_CURRENT_BINARY_DIR}")

add_library(nacs-kernel SHARED
  ${nacs_kernel_SRCS})

target_link_libraries(nacs-kernel ${nacs_kernel_LINKS})

set_target_properties(nacs-kernel PROPERTIES
  VERSION 0.1
  SOVERSION 0
  COMPILE_FLAGS "-fvisibility=hidden"
  LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")

install(TARGETS nacs-kernel LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}")
install(FILES ${nacs_kernel_HDRS}
  DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/nacs-kernel")
