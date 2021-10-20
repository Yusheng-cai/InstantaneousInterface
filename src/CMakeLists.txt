add_library(IIsrc SHARED)

target_sources(IIsrc PRIVATE
    "AtomGroup.cpp"
    "AtomGroupParsingStrategy.cpp"
    "BoundingBox.cpp"
    "DensityField.cpp"
    "Driver.cpp"
    "Field.cpp"
    "SimulationBox.cpp"
    "SimulationState.cpp"
    "MarchingCubesWrapper.cpp"
    "AverageField.cpp"
    "LinAlgTools.cpp"
    "Mesh.cpp"
    "Curvature.cpp"
    "CurvatureFDM.cpp"
    "CurvatureTensor.cpp"
    "Graph.cpp"
    "CurvatureCurveFit.cpp"
    "MeshRefineStrategy.cpp"
    "MeshLoopRefine.cpp"
    "CurvatureJetFit.cpp"
)

set(IIsrc_HEADER
    "AtomGroup.h"
    "AtomGroupParsingStrategy.h"
    "BoundingBox.h"
    "DensityField.h"
    "Driver.h"
    "Field.h"
    "SimulationBox.h"
    "SimulationState.h"
    "GaussianCoarseGrainFunction.h"
    "MarchingCubesWrapper.h"
    "LinAlgTools.h"
    "AverageField.h"
    "Mesh.h"
    "Curvature.h"
    "CurvatureFDM.h"
    "CurvatureTensor.h"
    "Graph.h"
    "CurvatureCurveFit.h"
    "MeshRefineStrategy.h"
    "MeshLoopRefine.h"
    "CurvatureJetFit.h"
)

target_include_directories(IIsrc
                        PUBLIC 
                        $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}>
                        $<INSTALL_INTERFACE:${PROJECT_BINARY_DIR}/include>
)

target_link_libraries(IIsrc PUBLIC IIxdr IItools IIparallel MarchingCubes)

target_include_directories(IIsrc PUBLIC "/usr/local/include/")

include(InstallLibraryWithStandardSetup)
InstallLibraryWithStandardSetup(IIsrc 
    INCLUDE_SUBDIR IIsrc
    INCLUDES ${IIsrc_HEADER})