add_library(IIHLBFGS SHARED)

target_sources(IIHLBFGS PRIVATE
    "HLBFGS_BLAS.cpp"
    "HLBFGS.cpp"
    "ICFS.cpp"
    "LineSearch.cpp"
    "Lite_Sparse_Matrix.cpp"
)

set(IIHLBFGS_HEADER
    "HLBFGS_BLAS.h"
    "HLBFGS.h"
    "ICFS.h"
    "LineSearch.h"
    "Lite_Sparse_Matrix.h"
    "Sparse_Entry.h"
)

target_include_directories(IIHLBFGS
                        PUBLIC 
                        $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}>
                        $<INSTALL_INTERFACE:${PROJECT_BINARY_DIR}/include>
)

include(InstallLibraryWithStandardSetup)
InstallLibraryWithStandardSetup(IIHLBFGS
    INCLUDE_SUBDIR IIHLBFGS
    INCLUDES ${IIHLBFGS_HEADER})