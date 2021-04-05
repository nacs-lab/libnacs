cmake_minimum_required(VERSION 3.7)

if(POLICY CMP0079)
  # Causes yaml-cpp linking issue on msys2 otherwise.
  cmake_policy(SET CMP0079 NEW)
endif()

set(MAJOR_VERSION 10)
set(MINOR_VERSION 2)
set(PATCH_VERSION 0)

project(libnacs VERSION "${MAJOR_VERSION}.${MINOR_VERSION}.${PATCH_VERSION}"
  LANGUAGES C CXX ASM)

## Setup
set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH}"
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

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
option(ENABLE_INSTALL_NONLIB
  "Enable installation of files other than the libraries (e.g. headers and pkgconfig)"
  ${IS_UNIX})
option(ENABLE_SIMD "Enable simd optimization" On)
option(SEQ_IN_UTILS "Build nacs-seq as part of nacs-utils" Off)
option(ENABLE_TESTING "Enable building and running testing" On)
if (ENABLE_LLVM AND STATIC_LLVM)
  set(SEQ_IN_UTILS On FORCE)
endif()

include(CMakeCMacros)
include(GNUInstallDirs)
include(CheckCXXSourceCompiles)

set(INSTALL_PKGCONFIG_DIR "${CMAKE_INSTALL_FULL_LIBDIR}/pkgconfig"
  CACHE PATH "Installation directory for pkgconfig (.pc) files")

# Remove rdynamic
set(CMAKE_SHARED_LIBRARY_LINK_C_FLAGS)
set(CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED On)
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

find_package(PkgConfig REQUIRED)

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
  # Require this for now...
  find_package(yaml-cpp REQUIRED CONFIG NO_SYSTEM_ENVIRONMENT_PATH)
endif()

if(ENABLE_SIMD)
  find_package(Sleef REQUIRED)
  set(CMAKE_REQUIRED_INCLUDES ${SLEEF_INCLUDE_DIR})
  set(CMAKE_REQUIRED_LIBRARIES ${SLEEF_LIBRARIES})
  check_cxx_source_compiles("
      #include <sleef.h>
      int main()
      {
          ::Sleef_remainder;
      }" SLEEF_HAS_REMAINDER)
  if(SLEEF_HAS_REMAINDER)
    # For old Sleef version
    add_definitions(-DNACS_SLEEF_HAS_REMAINDER)
  endif()
endif()
find_package(OpenLibm)
pkg_check_modules(ZMQ REQUIRED IMPORTED_TARGET libzmq>=4.0)
if(APPLE)
  find_package(cppzmq)
  target_link_libraries(PkgConfig::ZMQ INTERFACE cppzmq)
endif()

add_definitions(-Wall -Wextra -g -pthread -fvisibility=hidden -Wconversion
  -fno-math-errno -Wno-sign-conversion -Wsign-compare
  -Wno-unused-parameter -Wno-psabi -D_FILE_OFFSET_BITS=64)
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  # We have GCC specific warning options that makes Clang unhappy.
  # Adding this unconditionally cause GCC to issue a warning if there are other warnings
  # so only add this conditionally for clang.
  add_definitions(-Wno-unknown-warning-option)
endif()

set(LINKER_FLAGS -pthread -fno-math-errno)
if(NOT APPLE)
  set(LINKER_FLAGS -Wl,--as-needed -Wl,--no-undefined -Wl,--gc-sections ${LINKER_FLAGS})
endif()

if(IS_WINDOWS)
  # GCC is somehow not happy about the `<inttypes.h>` we find even though it appears to work
  cmake_c_add_flags(CMAKE_C_FLAGS -Wno-format -Wno-format-extra-args)
  cmake_c_add_flags(CMAKE_CXX_FLAGS -Wno-format -Wno-format-extra-args)
  set(LINKER_FLAGS ${LINKER_FLAGS} -static-libstdc++)
  # The "runtime" part of all the libraries and executables
  # needs to be in the same directory on windows for them to be runnable
  # in the build directory due to the lack of rpath.
  # We don't need to set this for library and archive output since
  # those are only used during compile time and aren't subjected to
  # the windows runtime loader limitations.
  # Additionally, targets with different name prefix and suffix
  # (e.g. libxxx.dll and xxx.exe) may generate library output of the same name (libxxx.dll.a)
  # and setting a global library output directory may cause confilcts between them.
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/output-bin")
  # Workaround CMake bug for Clang
  set(CMAKE_ASM_COMPILE_OPTIONS_PIC)
  set(CMAKE_ASM_COMPILE_OPTIONS_PIE)
elseif(NOT APPLE)
  # Do no use `-Bsymbolic`.
  # The difference between that and `-Bsymbolic-functions` is that `-Bsymbolic`
  # also cause all the global variables to refer to the local version.
  # This is actually problematic when the global variable is accessed externally
  # since GCC can actually emit a static copy of the variable in a program
  # linked to the library and rely on the runtime linker
  # to overwrite the address in the library with the external one.
  # With `-Bsymbolic` the two will be refering to
  # two different versions of the global variable.
  set(LINKER_FLAGS ${LINKER_FLAGS} -Wl,-Bsymbolic-functions)
endif()

cmake_c_add_flags(CMAKE_C_FLAGS -std=gnu11)
cmake_c_add_flags(CMAKE_CXX_FLAGS -D__STDCPP_WANT_MATH_SPEC_FUNCS__)
cmake_c_add_flags(CMAKE_EXE_LINKER_FLAGS ${LINKER_FLAGS})
cmake_c_add_flags(CMAKE_STATIC_LINKER_FLAGS ${LINKER_FLAGS})
cmake_c_add_flags(CMAKE_SHARED_LINKER_FLAGS ${LINKER_FLAGS})
cmake_c_add_flags(CMAKE_MODULE_LINKER_FLAGS ${LINKER_FLAGS})
cmake_c_add_flags(CMAKE_C_FLAGS_DEBUG -ggdb3 -O0)
cmake_c_add_flags(CMAKE_CXX_FLAGS_DEBUG -ggdb3 -O0 -D_GLIBCXX_ASSERTIONS)
cmake_c_add_flags(CMAKE_C_FLAGS_RELEASE -O3)
cmake_c_add_flags(CMAKE_CXX_FLAGS_RELEASE -O3)
# GCC does not turn on loop vectorizer at -O2
cmake_c_add_flags(CMAKE_C_FLAGS_RELWITHDEBINFO -O3)
cmake_c_add_flags(CMAKE_CXX_FLAGS_RELWITHDEBINFO -O3)

include_directories("${CMAKE_CURRENT_SOURCE_DIR}/lib")
include_directories("${CMAKE_CURRENT_BINARY_DIR}/lib/nacs-utils") # config.h

if(IS_WINDOWS)
  # Used by test
  cmake_utils_cygpath_m(INSTALL_WINDOWS_BINDIR "${CMAKE_INSTALL_FULL_BINDIR}")
else()
  set(INSTALL_WINDOWS_BINDIR "${CMAKE_INSTALL_FULL_BINDIR}")
endif()

add_subdirectory(lib)
add_subdirectory(tools)
if(ENABLE_TESTING)
  enable_testing()
  add_subdirectory(test)
endif()

function(get_soname out prefix libs_var)
  set(${out} "" PARENT_SCOPE)
  if(CMAKE_OBJDUMP)
    execute_process(COMMAND "${CMAKE_OBJDUMP}" -p "${${libs_var}}"
      RESULT_VARIABLE objdump_res
      OUTPUT_VARIABLE objdump_out
      ERROR_QUIET)
    if(objdump_res EQUAL 0)
      string(REGEX MATCH "SONAME[^\n]*${prefix}\\.so\\.([0-9]*)" match "${objdump_out}")
      set(${out} "${CMAKE_MATCH_1}" PARENT_SCOPE)
    endif()
  endif()
endfunction()

if(EXISTS "${CMAKE_ROOT}/Modules/CPack.cmake")
  include(InstallRequiredSystemLibraries)

  set(CPACK_SET_DESTDIR "on")
  set(CPACK_PACKAGING_INSTALL_PREFIX "/tmp")
  set(CPACK_GENERATOR "DEB")

  set(CPACK_PACKAGE_NAME "${CMAKE_PROJECT_NAME}" CACHE STRING "DEB package name")
  set(CPACK_PACKAGE_DESCRIPTION "Libraries for NaCs lab control")
  set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Libraries for NaCs lab control")
  set(CPACK_PACKAGE_VENDOR "Harvard")
  set(CPACK_PACKAGE_CONTACT "Yichao Yu")
  if("${CPACK_DEBIAN_PACKAGE_RELEASE}" STREQUAL "")
    set(DEB_SUFFIX "")
  else()
    set(DEB_SUFFIX "-${CPACK_DEBIAN_PACKAGE_RELEASE}")
  endif()
  string(TOLOWER "${CMAKE_SYSTEM_PROCESSOR}" system_processor)
  if("${system_processor}" STREQUAL "x86_64" OR
      "${system_processor}" STREQUAL "x86-64" OR
      "${system_processor}" STREQUAL "amd64")
    set(CPACK_DEBIAN_PACKAGE_ARCHITECTURE "amd64")
  elseif("${system_processor}" STREQUAL "x86" OR
      "${system_processor}" STREQUAL "i386" OR
      "${system_processor}" STREQUAL "i486" OR
      "${system_processor}" STREQUAL "i586" OR
      "${system_processor}" STREQUAL "i686")
    set(CPACK_DEBIAN_PACKAGE_ARCHITECTURE "i386")
  elseif("${system_processor}" STREQUAL "arm" OR
      "${system_processor}" STREQUAL "armv7l" OR
      "${system_processor}" STREQUAL "armhf")
    set(CPACK_DEBIAN_PACKAGE_ARCHITECTURE "armhf")
  elseif("${system_processor}" STREQUAL "aarch64" OR
      "${system_processor}" STREQUAL "arm64")
    set(CPACK_DEBIAN_PACKAGE_ARCHITECTURE "arm64")
  else()
    message(WARNING "Unrecognized architecture ${system_processor}")
    set(CPACK_DEBIAN_PACKAGE_ARCHITECTURE "${system_processor}")
  endif()
  set(CPACK_PACKAGE_FILE_NAME
    "${CPACK_PACKAGE_NAME}_${MAJOR_VERSION}.${MINOR_VERSION}.${PATCH_VERSION}${DEB_SUFFIX}_${CPACK_DEBIAN_PACKAGE_ARCHITECTURE}")

  set(CPACK_DEBIAN_PACKAGE_DEPENDS "libstdc++6 (>= 4.2.1)")
  get_soname(OPENLIBM_SOVER libopenlibm OPENLIBM_LIBRARIES)
  set(CPACK_DEBIAN_PACKAGE_DEPENDS
    "${CPACK_DEBIAN_PACKAGE_DEPENDS}, libopenlibm${OPENLIBM_SOVER}")
  if(ENABLE_LLVM AND NOT STATIC_LLVM)
    string(REPLACE "." ";" LLVM_PACKAGE_VERSION_LIST ${LLVM_PACKAGE_VERSION})
    list(GET LLVM_PACKAGE_VERSION_LIST 0 LLVM_PACKAGE_VERSION_MAJOR)
    set(CPACK_DEBIAN_PACKAGE_DEPENDS
      "${CPACK_DEBIAN_PACKAGE_DEPENDS}, libllvm${LLVM_PACKAGE_VERSION_MAJOR}")

    string(REGEX MATCH "^[0-9]*\\.[0-9]*" YAML_CPP_PKGVER "${yaml-cpp_VERSION}")
    if("${YAML_CPP_PKGVER}" STREQUAL "0.3" OR "${YAML_CPP_PKGVER}" STREQUAL "0.5")
      # This is the suffix in the Ubuntu 18.04 and Debian stretch package.
      # I'm not really sure why it has a `v5` suffix but hopefully
      # this is good enough for the ubuntu versions we support.
      # If this is not good enough, we might need a separate option
      # to pass the correct values in...
      set(YAML_CPP_PKGVER "${YAML_CPP_PKGVER}v5")
    endif()
    set(CPACK_DEBIAN_PACKAGE_DEPENDS
      "${CPACK_DEBIAN_PACKAGE_DEPENDS}, libyaml-cpp${YAML_CPP_PKGVER}")
  endif()
  if(ENABLE_SIMD)
    get_soname(SLEEF_SOVER libsleef SLEEF_LIBRARIES)
    set(CPACK_DEBIAN_PACKAGE_DEPENDS
      "${CPACK_DEBIAN_PACKAGE_DEPENDS}, libsleef${SLEEF_SOVER}")
  endif()
  # FIXME: The ZMQ version detection doesn't seem to work reliably on ubuntu
  # so let's just hard code it for now...
  # get_soname(ZMQ_SOVER libzmq ZMQ_LIBRARIES)
  set(ZMQ_SOVER 5)
  set(CPACK_DEBIAN_PACKAGE_DEPENDS
    "${CPACK_DEBIAN_PACKAGE_DEPENDS}, libzmq${ZMQ_SOVER}")

  set(CPACK_DEBIAN_PACKAGE_PRIORITY "optional")
  set(CPACK_DEBIAN_PACKAGE_SECTION "")

  set(CPACK_COMPONENTS_ALL Libraries ApplicationData)
  include(CPack)
endif()
