#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "gamera::gamera" for configuration "Release"
set_property(TARGET gamera::gamera APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(gamera::gamera PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libgamera.so"
  IMPORTED_SONAME_RELEASE "libgamera.so"
  )

list(APPEND _cmake_import_check_targets gamera::gamera )
list(APPEND _cmake_import_check_files_for_gamera::gamera "${_IMPORT_PREFIX}/lib/libgamera.so" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
