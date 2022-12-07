add_library(IIICP SHARED)

target_sources(IIICP PRIVATE
    "ICP.cpp"
    "svd.cpp"
    "KdTree.cpp"
)

set(IIsrc_HEADER
    "ICP.h"
    "KdTree.h"
    "defs_and_types.h"
    "svd.h"
)

target_include_directories(IIICP
                        PUBLIC 
                        $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}>
                        $<INSTALL_INTERFACE:${PROJECT_BINARY_DIR}/include>
)

include(InstallLibraryWithStandardSetup)
InstallLibraryWithStandardSetup(IIICP
    INCLUDE_SUBDIR IIICP
    INCLUDES ${IIICP_HEADER})
