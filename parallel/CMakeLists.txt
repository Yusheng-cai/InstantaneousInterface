add_library(IIparallel SHARED)

target_sources(IIparallel PRIVATE
    "OpenMP.cpp"
)

set(IIPARALLEL_HEADER
    "MPI_environment.h"
    "OpenMP_buffer.h"
    "OpenMP.h")

target_include_directories(IIparallel PUBLIC 
                         "${CMAKE_SOURCE_DIR}")

# target_compile_options(parallel PRIVATE -fPIC)

message("Required package of OpenMP is ${OpenMP_CXX_FOUND}")

# resource: https://gitlab.kitware.com/cmake/cmake/-/issues/21818
if(OpenMP_CXX_FOUND)
    target_link_libraries(IIparallel PUBLIC OpenMP::OpenMP_CXX)
endif()

include(InstallLibraryWithStandardSetup)
InstallLibraryWithStandardSetup(IIparallel 
    INCLUDE_SUBDIR IIparallel 
    INCLUDES ${IIPARALLEL_HEADER})
    