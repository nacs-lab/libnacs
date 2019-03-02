#!/bin/sh

cmake="$1"
get_soname="$2"
lib_suffix="$3"
out_file="$4"

write_file() {
    printf "#define %s \"%s\"\\n" "$1" "$("$get_soname" "$2")"
}

{
    write_file NACS_OPENLIBM_NAME "libopenlibm$lib_suffix"
} > "${out_file}.tmp"

"$cmake" -E copy_if_different "${out_file}.tmp" "${out_file}"
