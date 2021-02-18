add_library(xrd STATIC
    ${CMAKE_CURRENT_SOURCE_DIR}/basis.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/crystal.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/diffraction.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/lattice.cpp)
set_target_properties(xrd PROPERTIES
    CXX_VISIBILITY_PRESET "hidden")
target_link_libraries(xrd
    PUBLIC
    xrd_utils
    PRIVATE
    OpenMP::OpenMP_CXX)

add_executable(xrd_simulation
    ${CMAKE_CURRENT_SOURCE_DIR}/main_simulation.cpp)
set_target_properties(xrd_simulation PROPERTIES
    CXX_VISIBILITY_PRESET "hidden")
target_link_libraries(xrd_simulation
    PRIVATE
    OpenMP::OpenMP_CXX
    xrd)

add_executable(xrd_optimisation
    ${CMAKE_CURRENT_SOURCE_DIR}/main_optimisation.cpp)
set_target_properties(xrd_optimisation PROPERTIES
    CXX_VISIBILITY_PRESET "hidden")
target_link_libraries(xrd_optimisation
    PRIVATE
    xrd)

add_executable(xrd_background
    ${CMAKE_CURRENT_SOURCE_DIR}/main_background.cpp)
set_target_properties(xrd_background PROPERTIES
    CXX_VISIBILITY_PRESET "hidden")
target_link_libraries(xrd_background
    PRIVATE
    xrd_utils)

add_executable(xrd_peaks
    ${CMAKE_CURRENT_SOURCE_DIR}/main_peak_finder.cpp)
  set_target_properties(xrd_peaks PROPERTIES
      CXX_VISIBILITY_PRESET "hidden")
  target_link_libraries(xrd_peaks
      PRIVATE
        xrd_utils)