project(LibNaCs)

## Setup
cmake_minimum_required(VERSION 2.8)
set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH}"
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

include(CMakeCMacros)
include(GNUInstallDirs)

# Remove rdynamic
set(CMAKE_SHARED_LIBRARY_LINK_C_FLAGS)
set(CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS)

add_definitions(-Wall -Wextra -g -pthread -fvisibility=hidden
  -O2 -fno-math-errno -Wconversion -Wno-sign-conversion)

if(NOT CMAKE_C_COMPILER MATCHES clang)
  add_definitions(-flto)
endif()

cmake_c_fix_include_path(lib/utils nacs-utils)
cmake_c_fix_include_path(lib/kernel nacs-kernel)

set(LINKER_FLAGS -Wl,--as-needed -Wl,--no-undefined -Wl,--gc-sections
  -pthread -O2 -fno-math-errno)

if(NOT CMAKE_C_COMPILER MATCHES clang)
  set(LINKER_FLAGS ${LINKER_FLAGS} -flto)
endif()

cmake_c_add_flags(CMAKE_C_FLAGS -std=gnu11)
cmake_c_add_flags(CMAKE_CXX_FLAGS -std=gnu++14)
cmake_c_add_flags(CMAKE_EXE_LINKER_FLAGS ${LINKER_FLAGS})
cmake_c_add_flags(CMAKE_STATIC_LINKER_FLAGS ${LINKER_FLAGS})
cmake_c_add_flags(CMAKE_SHARED_LINKER_FLAGS ${LINKER_FLAGS})
cmake_c_add_flags(CMAKE_MODULE_LINKER_FLAGS ${LINKER_FLAGS})

cmake_c_include_fix_path()

add_subdirectory(lib)
add_subdirectory(test)
