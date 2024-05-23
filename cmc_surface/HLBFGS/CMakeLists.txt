add_library(HLBFGS SHARED)

target_sources(HLBFGS PRIVATE
    "HLBFGS_BLAS.cpp"
    "HLBFGS.cpp"
    "ICFS.cpp"
    "LineSearch.cpp"
    "Lite_Sparse_Matrix.cpp"
)

set(HLBFGS_HEADER
    "HLBFGS_BLAS.h"
    "HLBFGS.h"
    "ICFS.h"
    "LineSearch.h"
    "Lite_Sparse_Matrix.h"
    "Sparse_Entry.h"
)

target_include_directories(HLBFGS
                        PUBLIC 
                        $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}>
                        $<INSTALL_INTERFACE:${PROJECT_BINARY_DIR}/include>
)