add_library(IItools SHARED)

target_sources(IItools PRIVATE
    "RandomNumber.cpp"
    "InputParser.cpp"
    "CommandLineArguments.cpp"
    "FileSystem.cpp"
    "OutputFunction.cpp"
)

set(IITOOLS_HEADER
    "Assert.h"
    "CommonTypes.h"
    "GenericFactory.h"
    "InputParser.h"
    "RandomNumber.h"
    "CommandLineArguments.h"
    "Registry.h"
    "FileSystem.h"
    "OutputFunction.h"
)

target_include_directories(IItools
                        PUBLIC 
                        $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}>
                        $<INSTALL_INTERFACE:${PROJECT_BINARY_DIR}/include>
)
# target_compile_options(tools PRIVATE -fPIC)


include(InstallLibraryWithStandardSetup)
InstallLibraryWithStandardSetup(IItools 
    INCLUDE_SUBDIR IItools
    INCLUDES ${IITOOLS_HEADER})

