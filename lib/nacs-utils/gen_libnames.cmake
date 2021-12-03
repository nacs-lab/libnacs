#

# input variables
# CMAKE_OBJDUMP
# OUTPUT_FILE
# NAMES
# <NAMES>

function(write_barename name file)
  get_filename_component(libname "${file}" NAME)
  file(APPEND "${OUTPUT_FILE}.tmp" "#define NACS_${name}_NAME \"${libname}\"\n")
endfunction()

function(write_file name file)
  if(APPLE)
    # As a hack, we include the full path on macOS since I haven't yet found
    # a better solution to make finding homebrew libraries easier...
    file(APPEND "${OUTPUT_FILE}.tmp" "#define NACS_${name}_NAME \"${file}\"\n")
    return()
  endif()
  if(NOT CMAKE_OBJDUMP)
    write_barename("${name}" "${file}")
    return()
  endif()
  execute_process(COMMAND "${CMAKE_OBJDUMP}" -p "${file}"
    RESULT_VARIABLE objdump_res
    OUTPUT_VARIABLE objdump_out
    ERROR_QUIET)
  if(NOT objdump_res EQUAL 0)
    write_barename("${name}" "${file}")
    return()
  endif()
  string(REGEX MATCH "SONAME[ \t\n]*([^ \t\n]*)" match "${objdump_out}")
  write_barename("${name}" "${CMAKE_MATCH_1}")
endfunction()

file(WRITE "${OUTPUT_FILE}.tmp" "")

foreach(name ${NAMES})
  write_file(${name} "${${name}}")
endforeach()

# This is copy_if_different
configure_file("${OUTPUT_FILE}.tmp" "${OUTPUT_FILE}" COPYONLY)
