add_library(IIModelViewer SHARED)

target_sources(IIModelViewer PRIVATE
    "Shader.cpp"
    "glad.c"
)

set(IIModelViewer_HEADER
    "Shader.h"
    "glad.h"
)

target_include_directories(IIModelViewer
                        PUBLIC 
                        $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}>
                        $<INSTALL_INTERFACE:${PROJECT_BINARY_DIR}/include>
)

include(InstallLibraryWithStandardSetup)
InstallLibraryWithStandardSetup(IIModelViewer
    INCLUDE_SUBDIR IIModelViewer
    INCLUDES ${IIModelViewer_HEADER})
