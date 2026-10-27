cmake_minimum_required(VERSION 3.17...3.24 FATAL_ERROR)

project(
    GAMERA VERSION 1.0.0 LANGUAGES CXX
)

include(GNUInstallDirs)

option(BUILD_SHARED_LIBS "Build shared libraries" ON)
find_package(GSL REQUIRED)

add_library(gamera)
add_library(gamera::gamera ALIAS gamera)

target_link_libraries(gamera PUBLIC GSL::gsl)
target_include_directories(gamera PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

# c++11 without extensions
target_compile_features(gamera PUBLIC cxx_std_11)
set_target_properties(gamera PROPERTIES CXX_EXTENSIONS OFF CXX_STANDARD_REQUIRED ON)

add_subdirectory(src)

# installation / exports
install(TARGETS gamera EXPORT ${PROJECT_NAME}Targets COMPONENT Runtime)
install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} COMPONENT Development)

# Creating the export files, so other projects can use this project via `find_package`
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
    ${PROJECT_NAME}ConfigVersion.cmake
    VERSION ${PACKAGE_VERSION}
    COMPATIBILITY SameMajorVersion
)

# create the file describing all targets exported by this project
install(
    EXPORT ${PROJECT_NAME}Targets
    FILE ${PROJECT_NAME}Targets.cmake
    NAMESPACE gamera::
    DESTINATION ${CMAKE_INSTALL_DATADIR}/cmake/${PROJECT_NAME}
)

# install the main Config file and the ConfigVersion file
install(
    FILES "${PROJECT_NAME}Config.cmake" "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
    DESTINATION ${CMAKE_INSTALL_DATADIR}/cmake/${PROJECT_NAME}
)
