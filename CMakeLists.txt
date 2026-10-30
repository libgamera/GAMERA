cmake_minimum_required(VERSION 3.10)
project(GAMERA VERSION 1.0 LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Set output directories
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

# Find GSL library
find_package(GSL REQUIRED)

# Include directories
include_directories(
    ${CMAKE_SOURCE_DIR}/include
    ${GSL_INCLUDE_DIRS}
)

# Source files
set(GAMERA_SOURCES
    src/Radiation.C
    src/Particles.C
    src/Utils.C
    src/Astro.C
    src/2D_interp/interp2d.C
    src/2D_interp/interp2d_spline.C
    src/2D_interp/bilinear.C
    src/2D_interp/bicubic.C
)

# Create shared library
add_library(gamera SHARED ${GAMERA_SOURCES})

# Add alias for users who want to use this library from build directory
add_library(GAMERA::gamera ALIAS gamera)

# Link GSL libraries
target_link_libraries(gamera ${GSL_LIBRARIES})

# Set compilation flags
target_compile_options(gamera PRIVATE
    -Wall
    -O2
    -g
    -fPIC
)

# Set include directories for the library itself (for building)
target_include_directories(gamera PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include/gamera>
)

# Installation rules
install(TARGETS gamera
    EXPORT GameraTargets
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    RUNTIME DESTINATION bin
    INCLUDES DESTINATION include
)

# Install headers to include/gamera/ (preserving directory structure)
install(DIRECTORY include/
    DESTINATION include/gamera
    FILES_MATCHING PATTERN "*.h"
)

# Create and install config files for easy consumption by other projects
install(EXPORT GameraTargets
    FILE GameraTargets.cmake
    NAMESPACE GAMERA::
    DESTINATION lib/cmake/GAMERA
)

# Create config file
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/GameraConfigVersion.cmake"
    VERSION ${GAMERA_VERSION}
    COMPATIBILITY AnyNewerVersion
)

install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/GameraConfigVersion.cmake"
    DESTINATION lib/cmake/GAMERA
)

# Create main config file
configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/GameraConfig.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/GameraConfig.cmake"
    @ONLY
)

install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/GameraConfig.cmake"
    DESTINATION lib/cmake/GAMERA
)

# Print configuration summary
message(STATUS "GAMERA Configuration:")
message(STATUS "  Version: ${GAMERA_VERSION}")
message(STATUS "  C++ Standard: ${CMAKE_CXX_STANDARD}")
message(STATUS "  Install prefix: ${CMAKE_INSTALL_PREFIX}")
message(STATUS "  GSL include dirs: ${GSL_INCLUDE_DIRS}")
message(STATUS "  GSL libraries: ${GSL_LIBRARIES}")