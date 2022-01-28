cmake_minimum_required(VERSION 3.18.1)

# Include the cmake scripts directory
LIST( APPEND CMAKE_MODULE_PATH
  ${CMAKE_SOURCE_DIR}/cmake/
)

# create the new project
project(InstantaneousInterface)

# Require C++11
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g -DMY_DEBUG")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -ffast-math")

# Requires OpenMP
find_package(OpenMP REQUIRED)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

#set( BOOST_ROOT "" CACHE PATH "Boost library path" )
#set( Boost_NO_SYSTEM_PATHS on CACHE BOOL "Do not search system for Boost" )
#find_package(Boost REQUIRED)
#message("Boost found = ${Boost_FOUND}")

find_package(CGAL)
if (CGAL_FOUND)
  message("CGAL FOUND.")
  message("CGAL INCLUDE DIR = ${CGAL_INCLUDE_DIRS}")
  message("CGAL LIBRARIES = ${CGAL_LIBRARIES}")
endif()
message("CMAKE ROOT = ${CMAKE_ROOT}")

# install the libraries as well as the executables in lib & bin folder in build/
# This has to be after project is defined
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/bin")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/lib")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/lib")

enable_testing()

add_subdirectory(src)
add_subdirectory(tools)
add_subdirectory(parallel)
add_subdirectory(xdr)
add_subdirectory(MarchingCubes)
add_subdirectory(test)

add_executable(InstantaneousInterface main.cpp)
add_executable(Curve curve.cpp)
add_executable(mesh_gen mesh_gen.cpp)

target_include_directories(InstantaneousInterface PUBLIC ${CMAKE_SOURCE_DIR})
target_include_directories(Curve PUBLIC ${CMAKE_SOURCE_DIR})
target_include_directories(mesh_gen PUBLIC ${CMAKE_SOURCE_DIR})

target_link_libraries(InstantaneousInterface PUBLIC IIxdr IIsrc IItools IIparallel MarchingCubes)
target_link_libraries(Curve PUBLIC IIxdr IIsrc IItools IIparallel MarchingCubes)
target_link_libraries(mesh_gen PUBLIC IIxdr IIsrc IItools IIparallel MarchingCubes)


install(TARGETS InstantaneousInterface DESTINATION bin) 
install(TARGETS Curve DESTINATION bin)
install(TARGETS mesh_gen DESTINATION bin)
