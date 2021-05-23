add_library(gslpp STATIC
    ${CMAKE_CURRENT_SOURCE_DIR}/error.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/integration.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/interp.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/spline.cpp)
  target_link_libraries(gslpp
      PRIVATE
        fmt::fmt
        PkgConfig::GSL)