cmake_minimum_required(VERSION 3.17...3.24 FATAL_ERROR)

project(
    GAMERA VERSION 1.0.0 LANGUAGES CXX
)


option(BUILD_SHARED_LIBS "Build shared libraries" ON)
option(PYTHON "Build gappa python bindings" ON)
find_package(GSL REQUIRED)

add_library(gamera)
add_library(gamera::gamera ALIAS gamera)
target_include_directories(gamera PUBLIC include)
target_link_libraries(gamera PUBLIC GSL::gsl)

# c++11 without extensions
target_compile_features(gamera PUBLIC cxx_std_11)
set_target_properties(gamera PROPERTIES CXX_EXTENSIONS OFF CXX_STANDARD_REQUIRED ON)

add_subdirectory(src)

install(TARGETS gamera)
