add_library(src SHARED)

target_sources(src PRIVATE
    "AtomGroup.cpp"
    "AtomGroupParsingStrategy.cpp"
    "BoundingBox.cpp"
    "DensityField.cpp"
    "Driver.cpp"
    "Field.cpp"
    "SimulationBox.cpp"
    "SimulationState.cpp"
)

set(src_HEADER
    "AtomGroup.h"
    "AtomGroupParsingStrategy.h"
    "BoundingBox.h"
    "DensityField.h"
    "Driver.h"
    "Field.h"
    "SimulationBox.h"
    "SimulationState.h"
    "GaussianCoarseGrainFunction.h"
)

target_include_directories(src
                        PUBLIC 
                        $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}>
                        $<INSTALL_INTERFACE:${PROJECT_BINARY_DIR}/include>
)

target_link_libraries(src PUBLIC xdr tools parallel)

include(InstallLibraryWithStandardSetup)
InstallLibraryWithStandardSetup(src 
    INCLUDE_SUBDIR src
    INCLUDES ${src_HEADER})