# Try to find Openlibm functionality
# Once done this will define
#
#  OPENLIBM_FOUND - system has Openlibm
#  OPENLIBM_LIBRARIES - Libraries needed to use Openlibm
#

find_library(OPENLIBM_LIBRARIES NAMES openlibm)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(OpenLibm DEFAULT_MSG OPENLIBM_LIBRARIES)

mark_as_advanced(OPENLIBM_INCLUDE_DIR OPENLIBM_LIBRARIES)
