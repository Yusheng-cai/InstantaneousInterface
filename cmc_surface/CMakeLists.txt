add_library(IIcmc_surface SHARED)

target_sources(IIcmc_surface PRIVATE
    "fast_rdt.cpp"
)

set(IICMC_SURFACE_HEADER
    "fast_rdt.h"
    "geex_utils.h"
    "edgemesh.h")

target_include_directories(IIcmc_surface PUBLIC 
                         "${CMAKE_SOURCE_DIR}")

add_subdirectory(HLBFGS)
target_link_libraries(IIcmc_surface PUBLIC IIHLBFGS)

include(InstallLibraryWithStandardSetup)
InstallLibraryWithStandardSetup(IIcmc_surface 
    INCLUDE_SUBDIR IIcmc_surface
    INCLUDES ${IICMC_SURFACE_HEADER})
    