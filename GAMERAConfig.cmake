# export file, so other projects can use this via `find_package`
include(CMakeFindDependencyMacro)

# we need to find the external dependencies here
# essentially we need to repeat all the `find_package` calls from the build
# here as `find_dependency`
find_dependency(GSL REQUIRED)


# Include the auto-generated targets file
include("${CMAKE_CURRENT_LIST_DIR}/GAMERATargets.cmake")
