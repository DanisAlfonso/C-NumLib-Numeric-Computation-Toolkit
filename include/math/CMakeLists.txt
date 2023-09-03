add_subdirectory(LinearAlgebra)
add_library(math_utils INTERFACE)
target_include_directories(math_utils INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
