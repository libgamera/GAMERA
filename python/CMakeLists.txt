cmake_minimum_required(VERSION 3.17...3.24 FATAL_ERROR)

project(GAPPA VERSION 1.0.0 LANGUAGES CXX)

include(GNUInstallDirs)
set(CMAKE_INSTALL_LIBDIR ".")

# RPATH setting so that the build libraries can load each other
if(APPLE)
    set(CMAKE_INSTALL_RPATH "@loader_path")
elseif(NOT WIN32)
    set(CMAKE_INSTALL_RPATH "$ORIGIN")
endif()


add_subdirectory(.. ${CMAKE_CURRENT_BINARY_DIR}/gamera)

find_package(Python COMPONENTS Interpreter Development.Module REQUIRED)
find_package(SWIG COMPONENTS python REQUIRED)

include(UseSWIG)

set_property(SOURCE gappa.i PROPERTY CPLUSPLUS ON)
swig_add_library(
    gappa
    LANGUAGE python
    OUTPUT_DIR "${SKBUILD_PLATLIB_DIR}"
    SOURCES gappa.i
)

target_link_libraries(gappa PRIVATE gamera::gamera)

if(WIN32)
  set_property(TARGET gappa PROPERTY SUFFIX ".${Python_SOABI}.pyd")
else()
  set_property(TARGET gappa
               PROPERTY SUFFIX ".${Python_SOABI}${CMAKE_SHARED_MODULE_SUFFIX}")
endif()
target_link_libraries(gappa PRIVATE Python::Module)

install(TARGETS gappa DESTINATION . COMPONENT Python)
