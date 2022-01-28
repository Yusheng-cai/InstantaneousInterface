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
    "UmbrellaSmoothing.cpp"
    "MeshCurvatureflow.cpp"
    "MeshCurvature.cpp"
    "Registry.cpp"
    "KDTree.cpp"
    "marching_cubes.cpp"
    "MeshGeneration.cpp"
    "Mesh2d.cpp"
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
    "UmbrellaSmoothing.h"
    "MeshCurvatureflow.h"
    "MeshCurvature.h"
    "Registry.h"
    "KDTree.h"
    "marching_cubes.hpp"
    "MeshGeneration.h"
    "Mesh2d.h"
)

target_include_directories(IIsrc
                        PUBLIC 
                        $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}>
                        $<INSTALL_INTERFACE:${PROJECT_BINARY_DIR}/include>
)


#if(Boost_FOUND)
#    target_include_directories(IIsrc PUBLIC ${Boost_INCLUDE_DIRS})
#endif()

if (CGAL_FOUND)
   target_include_directories(IIsrc PUBLIC ${CGAL_INCLUDE_DIRS})
endif()
target_link_libraries(IIsrc PUBLIC IIxdr IItools IIparallel MarchingCubes CGAL::CGAL)
#message("Boost include dir = ${Boost_INCLUDE_DIRS}")
#message("Boost library dir = ${Boost_LIBRARIES}")

include(InstallLibraryWithStandardSetup)
InstallLibraryWithStandardSetup(IIsrc 
    INCLUDE_SUBDIR IIsrc
    INCLUDES ${IIsrc_HEADER})