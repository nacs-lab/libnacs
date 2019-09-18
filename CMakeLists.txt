cmake_minimum_required(VERSION 3.0)

project(libnacs C CXX ASM)

## Setup
set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH}"
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

set(MAJOR_VERSION 9)
set(MINOR_VERSION 0)
set(PATCH_VERSION 0)

if(WIN32 OR CYGWIN OR MINGW OR MSYS)
  set(IS_WINDOWS True)
  set(IS_UNIX False)
else()
  set(IS_WINDOWS False)
  set(IS_UNIX True)
endif()

# Options
option(ENABLE_KNACS "Enable kernel driver" Off)
option(ENABLE_LLVM "Enable LLVM JIT" Off)
option(STATIC_LLVM "Statically link LLVM" Off)
option(ENABLE_PYWRAPPER "Enable python wrapper for nacs-seq" On)
option(ENABLE_INSTALL_NONLIB
  "Enable installation of files other than the libraries (e.g. headers and pkgconfig)"
  ${IS_UNIX})
option(ENABLE_SIMD "Enable simd optimization" On)

include(CMakeCMacros)
include(GNUInstallDirs)
include(CheckCXXSourceCompiles)

set(INSTALL_PKGCONFIG_DIR "${CMAKE_INSTALL_FULL_LIBDIR}/pkgconfig"
  CACHE PATH "Installation directory for pkgconfig (.pc) files")

# Remove rdynamic
set(CMAKE_SHARED_LIBRARY_LINK_C_FLAGS)
set(CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS)

set(CMAKE_REQUIRED_FLAGS "-std=gnu++14")
check_cxx_source_compiles("
      #include <math.h>
      int main()
      {
          ::exp10;
      }" HAS_EXP10)
if(HAS_EXP10)
  # For the interpreter on windows
  add_definitions(-DNACS_HAS_EXP10)
endif()

if(ENABLE_LLVM)
  # LLVM uses the path of the cmake file to determine the llvm install prefix.
  # Since `PATH` might contain `/bin` which will be reduced to `/`,
  # and since `/lib` is often symlinked to `/usr/lib` these days,
  # cmake could find the llvm config file from `/lib` first instead of `/usr/lib`.
  # This causes the include path to be wrong (`/include` instead of `/usr/include`)
  # Disable searching though `PATH` to fix this.
  find_package(LLVM REQUIRED CONFIG NO_SYSTEM_ENVIRONMENT_PATH)

  message(STATUS "Found LLVM ${LLVM_PACKAGE_VERSION}")
  message(STATUS "Using LLVMConfig.cmake in: ${LLVM_DIR}")
endif()

if(ENABLE_SIMD)
  find_package(Sleef REQUIRED)
endif()
find_package(OpenLibm)

add_definitions(-Wall -Wextra -g -pthread -fvisibility=hidden -Wconversion
  -fno-math-errno -Wno-sign-conversion -Wno-psabi -D_FILE_OFFSET_BITS=64)
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  # We have GCC specific warning options that makes Clang unhappy.
  # Adding this unconditionally cause GCC to issue a warning if there are other warnings
  # so only add this conditionally for clang.
  add_definitions(-Wno-unknown-warning-option)
endif()

cmake_c_fix_include_path(lib/utils nacs-utils)
cmake_c_fix_include_path(lib/kernel nacs-kernel)
cmake_c_fix_include_path(lib/seq nacs-seq)

set(LINKER_FLAGS -Wl,--as-needed -Wl,--no-undefined -Wl,--gc-sections -pthread -fno-math-errno)

if(IS_WINDOWS)
  # GCC is somehow not happy about the `<inttypes.h>` we find even though it appears to work
  cmake_c_add_flags(CMAKE_C_FLAGS -Wno-format -Wno-format-extra-args)
  cmake_c_add_flags(CMAKE_CXX_FLAGS -Wno-format -Wno-format-extra-args)
  set(LINKER_FLAGS ${LINKER_FLAGS} -static-libstdc++)
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/output-bin")
  set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/output-lib")
  set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/output-lib")
  # Workaround CMake bug for Clang
  set(CMAKE_ASM_COMPILE_OPTIONS_PIC)
  set(CMAKE_ASM_COMPILE_OPTIONS_PIE)
else()
  set(LINKER_FLAGS ${LINKER_FLAGS} -Wl,-Bsymbolic)
endif()

cmake_c_add_flags(CMAKE_C_FLAGS -std=gnu11)
cmake_c_add_flags(CMAKE_CXX_FLAGS -std=gnu++14 -fno-rtti -D__STDCPP_WANT_MATH_SPEC_FUNCS__)
cmake_c_add_flags(CMAKE_EXE_LINKER_FLAGS ${LINKER_FLAGS})
cmake_c_add_flags(CMAKE_STATIC_LINKER_FLAGS ${LINKER_FLAGS})
cmake_c_add_flags(CMAKE_SHARED_LINKER_FLAGS ${LINKER_FLAGS})
cmake_c_add_flags(CMAKE_MODULE_LINKER_FLAGS ${LINKER_FLAGS})
cmake_c_add_flags(CMAKE_C_FLAGS_DEBUG -ggdb3 -O0)
cmake_c_add_flags(CMAKE_CXX_FLAGS_DEBUG -ggdb3 -O0)

cmake_c_include_fix_path()

include_directories("${CMAKE_CURRENT_BINARY_DIR}/lib/utils") # config.h

if(IS_WINDOWS)
  # Used by python and test
  cmake_utils_cygpath_m(INSTALL_WINDOWS_BINDIR "${CMAKE_INSTALL_FULL_BINDIR}")
else()
  set(INSTALL_WINDOWS_BINDIR "${CMAKE_INSTALL_FULL_BINDIR}")
endif()

add_subdirectory(lib)
if(ENABLE_PYWRAPPER)
  add_subdirectory(python)
endif()
add_subdirectory(tools)
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

  set(CPACK_DEBIAN_PACKAGE_DEPENDS "libstdc++6 (>= 4.2.1), libopenlibm2")

  set(CPACK_DEBIAN_PACKAGE_PRIORITY "optional")
  set(CPACK_DEBIAN_PACKAGE_SECTION "")
  set(CPACK_DEBIAN_ARCHITECTURE ${CMAKE_SYSTEM_PROCESSOR})

  set(CPACK_COMPONENTS_ALL Libraries ApplicationData)
  include(CPack)
endif()
