add_library(IIModelViewer SHARED)

target_sources(IIModelViewer PRIVATE
    "Shader.cpp"
    "glad.c"
    "GLMesh.cpp"
)

set(IIModelViewer_HEADER
    "Shader.h"
    "glad.h"
    "GLMesh.h"
)

target_include_directories(IIModelViewer
                        PUBLIC 
                        $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}>
                        $<INSTALL_INTERFACE:${PROJECT_BINARY_DIR}/include>
)

target_link_libraries(IIsrc PUBLIC IIxdr IItools IIparallel CGAL::CGAL)

include(InstallLibraryWithStandardSetup)
InstallLibraryWithStandardSetup(IIModelViewer
    INCLUDE_SUBDIR IIModelViewer
    INCLUDES ${IIModelViewer_HEADER})
