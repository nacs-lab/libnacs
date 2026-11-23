add_definitions(-UNDEBUG)

if(IS_WINDOWS AND STATIC_LLVM)
  # Both Catch2 and the static LLVM libraries use `std::regex`.
  # The function-local static variables (and their guard variables) in it
  # may be put in sections with different names (`.data$` vs `.bss$`)
  # so the linker doesn't merge them and errors on multiple definitions.
  # These definitions are all identical so it's fine to pick any of them.
  cmake_c_add_flags(CMAKE_EXE_LINKER_FLAGS -Wl,--allow-multiple-definition)
endif()

find_package(Catch2 CONFIG)
if(NOT TARGET Catch2::Catch2)
  add_library(Catch2::Catch2 INTERFACE IMPORTED)
  set_target_properties(Catch2::Catch2 PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${CMAKE_CURRENT_SOURCE_DIR}/bundle")
endif()

# Compile the Catch2 implementation and `main` once and share it with all the tests.
add_library(nacs-catch-main STATIC catch_main.cpp)
target_link_libraries(nacs-catch-main PUBLIC Catch2::Catch2)
# Benchmarking changes the layout of Catch2's internal classes
# so it must be enabled consistently in all translation units.
target_compile_definitions(nacs-catch-main PUBLIC CATCH_CONFIG_ENABLE_BENCHMARKING)

configure_file(winpath_helper.h.in winpath_helper.h @ONLY)
include_directories("${CMAKE_CURRENT_BINARY_DIR}")

add_subdirectory(utils)
add_subdirectory(seq)
if(ENABLE_KNACS)
  add_subdirectory(kernel)
endif()
