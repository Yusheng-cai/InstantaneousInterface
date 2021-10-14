add_subdirectory(libxdr)
add_library(IIxdr SHARED)

target_sources(IIxdr PRIVATE
    "XdrWrapper.cpp"
    "TrrFile.cpp"
    "Frame.cpp"
    "XtcFile.cpp"
    "GroFile.cpp"
)

target_include_directories(IIxdr
                        PUBLIC ${CMAKE_SOURCE_DIR})


set(IIXDR_HEADER
    "Frame.h"
    "TrrFile.h"
    "XdrWrapper.h"
    "XtcFile.h"
    "GroFile.h"
)

target_link_libraries(IIxdr PRIVATE IItools xdrC)

include(InstallLibraryWithStandardSetup)
InstallLibraryWithStandardSetup(IIxdr 
    INCLUDE_SUBDIR IIxdr
    INCLUDES ${IIXDR_HEADER})