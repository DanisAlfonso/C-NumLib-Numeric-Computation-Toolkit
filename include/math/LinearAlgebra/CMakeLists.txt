add_library(LinearAlgebra INTERFACE)
target_sources(LinearAlgebra INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}/LUDecomposition.hpp
)
target_include_directories(LinearAlgebra INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(LinearAlgebra INTERFACE math_utils)
