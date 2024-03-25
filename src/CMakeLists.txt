add_library(IIsrc SHARED)

set(IIsrc_CPU_source 
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
    "CurvatureQuadraticFit.cpp"
    "Contour.cpp"
    "AverageFieldProperty.cpp"
    "AFP_shapes.cpp"
    "fsolve.cpp"
    "fsolve_wrapper.cpp"
)

set(IIsrc_GPU_source 
   "DensityFieldKernel.cu"
)

if (ENABLE_CUDA)
    message("Compiling CPU and CUDA.")
    target_sources(IIsrc PRIVATE
        ${IIsrc_CPU_source} ${IIsrc_GPU_source}
    )
else()
    message("Only compiling CPU code.")
    target_sources(IIsrc PRIVATE ${IIsrc_CPU_source})
endif()

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
    "CurvatureQuadraticFit.hpp"
    "Contour.h"
    "AverageFieldProperty.h"
    "AFP_shapes.h"
    "fsolve.hpp"
    "fsolve_wrapper.hpp"
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

if (ENABLE_CUDA)
   target_include_directories(IIsrc PUBLIC ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})
endif()

if (ENABLE_IGL)
    target_link_libraries(IIsrc PUBLIC IIxdr IItools IIparallel IIICP CGAL::CGAL igl::glfw)
else()
    target_link_libraries(IIsrc PUBLIC IIxdr IItools IIparallel IIICP CGAL::CGAL)
endif()


include(InstallLibraryWithStandardSetup)
InstallLibraryWithStandardSetup(IIsrc 
    INCLUDE_SUBDIR IIsrc
    INCLUDES ${IIsrc_HEADER})
