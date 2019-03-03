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

execute_process(COMMAND "${CMAKE_COMMAND}" -E copy_if_different
  "${OUTPUT_FILE}.tmp" "${OUTPUT_FILE}"
  OUTPUT_QUIET ERROR_QUIET)
