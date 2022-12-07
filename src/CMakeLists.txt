add_library(IIsrc SHARED)

target_sources(IIsrc PRIVATE
    "AtomGroup.cpp"
    "AtomGroupParsingStrategy.cpp"
    "Bin.cpp"
    "BoundingBox.cpp"
    "DensityField.cpp"
    "Driver.cpp"
    "Field.cpp"
    "SimulationBox.cpp"
    "SimulationState.cpp"
    "AverageField.cpp"
    "LinAlgTools.cpp"
    "Mesh.cpp"
    "Curvature.cpp"
    "CurvatureFDM.cpp"
    "CurvatureTensor.cpp"
    "Graph.cpp"
    "CurvatureCurveFit.cpp"
    "MeshRefineStrategy.cpp"
    "CurvatureJetFit.cpp"
    "UmbrellaSmoothing.cpp"
    "MeshCurvatureflow.cpp"
    "CurvatureEvolution.cpp"
    "Registry.cpp"
    "marching_cubes.cpp"
    "MeshGen2d.cpp"
    "MeshActions.cpp"
    "ShortEdgeRemoval.cpp"
    "ObtuseTriangleRemoval.cpp"
    "surface_nets.cpp"
    "LongEdgeRemoval.cpp"
)

set(IIsrc_HEADER
    "AtomGroup.h"
    "AtomGroupParsingStrategy.h"
    "Bin.h"
    "BoundingBox.h"
    "DensityField.h"
    "Driver.h"
    "Field.h"
    "SimulationBox.h"
    "SimulationState.h"
    "LinAlgTools.h"
    "AverageField.h"
    "Mesh.h"
    "Curvature.h"
    "CurvatureFDM.h"
    "CurvatureTensor.h"
    "Graph.h"
    "CurvatureCurveFit.h"
    "MeshRefineStrategy.h"
    "CurvatureJetFit.h"
    "UmbrellaSmoothing.h"
    "MeshCurvatureflow.h"
    "CurvatureEvolution.h"
    "Registry.h"
    "marching_cubes.hpp"
    "MeshGen2d.h"
    "MeshActions.h"
    "MeshPlaneIntersection.hpp"
    "ShortEdgeRemoval.h"
    "ObtuseTriangleRemoval.h"
    "surface_nets.h"
    "LongEdgeRemoval.h"
)

target_include_directories(IIsrc
                        PUBLIC 
                        $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}>
                        $<INSTALL_INTERFACE:${PROJECT_BINARY_DIR}/include>
)

add_subdirectory(ICP)

if (CGAL_FOUND)
   target_include_directories(IIsrc PUBLIC ${CGAL_INCLUDE_DIRS})
endif()
target_link_libraries(IIsrc PUBLIC IIxdr IItools IIparallel IIICP CGAL::CGAL)

include(InstallLibraryWithStandardSetup)
InstallLibraryWithStandardSetup(IIsrc 
    INCLUDE_SUBDIR IIsrc
    INCLUDES ${IIsrc_HEADER})
