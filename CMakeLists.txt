cmake_minimum_required(VERSION 3.17)
project(XRD)

set(CMAKE_CXX_STANDARD 20)

#find_package(cppcoro REQUIRED)
find_package (Eigen3 3.3 REQUIRED NO_MODULE)
find_package(fmt REQUIRED)
find_package(OpenMP REQUIRED COMPONENTS CXX)
find_package(range-v3 REQUIRED)

find_package(PkgConfig REQUIRED)
pkg_check_modules(GSL REQUIRED IMPORTED_TARGET gsl)

add_library(xrd_utils STATIC)
  set_target_properties(xrd_utils PROPERTIES
      CXX_VISIBILITY_PRESET         "hidden")
  target_compile_definitions(xrd_utils PUBLIC
      WITH_LOW_PRECISION)
  target_include_directories(xrd_utils PUBLIC
      ${CMAKE_CURRENT_SOURCE_DIR}/utils)
  target_sources(xrd_utils PRIVATE
      ${CMAKE_CURRENT_SOURCE_DIR}/utils/data/dataset_2d.cpp
      ${CMAKE_CURRENT_SOURCE_DIR}/utils/io.cpp
      ${CMAKE_CURRENT_SOURCE_DIR}/utils/math.cpp
      ${CMAKE_CURRENT_SOURCE_DIR}/utils/math/convolution.cpp
      ${CMAKE_CURRENT_SOURCE_DIR}/utils/math/peak_finder.cpp
      ${CMAKE_CURRENT_SOURCE_DIR}/utils/string.cpp
      ${CMAKE_CURRENT_SOURCE_DIR}/utils/timer.cpp)
  target_link_libraries(xrd_utils
      PUBLIC
        Eigen3::Eigen
        fmt::fmt
        gslpp
        nlohmann_json::nlohmann_json
        range-v3::range-v3
        PkgConfig::GSL
      PRIVATE
        OpenMP::OpenMP_CXX)

add_subdirectory(utils/gslpp)
add_subdirectory(xrd_simulation)

add_subdirectory(libs)