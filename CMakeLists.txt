project(libnacs)

## Setup
cmake_minimum_required(VERSION 2.8)
set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH}"
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

set(MAJOR_VERSION 2)
set(MINOR_VERSION 0)
set(PATCH_VERSION 0)

include(CMakeCMacros)
include(GNUInstallDirs)

set(INSTALL_PKGCONFIG_DIR "${CMAKE_INSTALL_FULL_LIBDIR}/pkgconfig"
  CACHE PATH "Installation directory for pkgconfig (.pc) files")

# Remove rdynamic
set(CMAKE_SHARED_LIBRARY_LINK_C_FLAGS)
set(CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS)

option(ENABLE_KNACS "Enable kernel driver" Off)

add_definitions(-Wall -Wextra -g -pthread -fvisibility=hidden
  -fno-math-errno -Wconversion -Wno-sign-conversion)

cmake_c_fix_include_path(lib/utils nacs-utils)
cmake_c_fix_include_path(lib/kernel nacs-kernel)
cmake_c_fix_include_path(lib/seq nacs-seq)

set(LINKER_FLAGS -Wl,--as-needed -Wl,--no-undefined -Wl,--gc-sections -pthread -fno-math-errno)

cmake_c_add_flags(CMAKE_C_FLAGS -std=gnu11)
cmake_c_add_flags(CMAKE_CXX_FLAGS -std=gnu++14)
cmake_c_add_flags(CMAKE_EXE_LINKER_FLAGS ${LINKER_FLAGS})
cmake_c_add_flags(CMAKE_STATIC_LINKER_FLAGS ${LINKER_FLAGS})
cmake_c_add_flags(CMAKE_SHARED_LINKER_FLAGS ${LINKER_FLAGS})
cmake_c_add_flags(CMAKE_MODULE_LINKER_FLAGS ${LINKER_FLAGS})

cmake_c_include_fix_path()

add_subdirectory(lib)
add_subdirectory(test)

if(EXISTS "${CMAKE_ROOT}/Modules/CPack.cmake")
  include(InstallRequiredSystemLibraries)

  set(CPACK_SET_DESTDIR "on")
  set(CPACK_PACKAGING_INSTALL_PREFIX "/tmp")
  set(CPACK_GENERATOR "DEB")

  set(CPACK_PACKAGE_DESCRIPTION "Libraries for NaCs lab control")
  set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Libraries for NaCs lab control")
  set(CPACK_PACKAGE_VENDOR "Harvard")
  set(CPACK_PACKAGE_CONTACT "Yichao Yu")
  set(CPACK_PACKAGE_VERSION_MAJOR "${MAJOR_VERSION}")
  set(CPACK_PACKAGE_VERSION_MINOR "${MINOR_VERSION}")
  set(CPACK_PACKAGE_VERSION_PATCH "${PATCH_VERSION}")
  set(CPACK_PACKAGE_FILE_NAME
    "${CMAKE_PROJECT_NAME}_${MAJOR_VERSION}.${MINOR_VERSION}.${CPACK_PACKAGE_VERSION_PATCH}")
  set(CPACK_SOURCE_PACKAGE_FILE_NAME
    "${CMAKE_PROJECT_NAME}_${MAJOR_VERSION}.${MINOR_VERSION}.${CPACK_PACKAGE_VERSION_PATCH}")

  set(CPACK_DEBIAN_PACKAGE_DEPENDS "libstdc++6 (>= 4.2.1)")

  set(CPACK_DEBIAN_PACKAGE_PRIORITY "optional")
  set(CPACK_DEBIAN_PACKAGE_SECTION "")
  set(CPACK_DEBIAN_ARCHITECTURE ${CMAKE_SYSTEM_PROCESSOR})

  set(CPACK_COMPONENTS_ALL Libraries ApplicationData)
  include(CPack)
endif()
