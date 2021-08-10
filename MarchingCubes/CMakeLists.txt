add_library(MarchingCubes SHARED)

target_sources(MarchingCubes PRIVATE
    "FindVertexes.cpp"
    "matrix.cpp"
    "PolygonizeCube.cpp"
    "TriangulateGlobal.cpp"
)

set(MarchingCubesHeader
    "matrix.hpp"
    "TriangulateGlobal.hpp"
    "TriangulateTab.hpp"
)

target_include_directories(MarchingCubes PUBLIC 
                         "${CMAKE_SOURCE_DIR}")

include(InstallLibraryWithStandardSetup)
InstallLibraryWithStandardSetup(parallel 
    INCLUDE_SUBDIR MarchingCubes 
    INCLUDES ${MarchingCubesHeader})
    