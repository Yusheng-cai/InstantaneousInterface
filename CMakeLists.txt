cmake_minimum_required(VERSION 3.18.1)

# Include the cmake scripts directory
LIST( PREPEND CMAKE_MODULE_PATH
  ${CMAKE_SOURCE_DIR}/cmake/
)

# create the new project
project(InstantaneousInterface LANGUAGES C CXX CUDA)

# libigl
include(libigl)
igl_include(glfw)

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

# set an option for using CUDA
option(ENABLE_CUDA "Using nvidia CUDA" OFF)

# Requires OpenGL
#find_package(OpenGL REQUIRED)
#if (OpenGL_FOUND)
#    message("OpenGL found")
#    message("lib = ${OPENGL_LIBRARIES}")
#    message("include = ${OPENGL_INCLUDE_DIR}")
#endif()

# Requires glut
#find_package(GLFW3 REQUIRED)
#if (GLFW3_FOUND)
#    message("GLFW FOUND")
#    message("GLFW Library directory = ${GLFW3_LIBRARY}")
#endif()

find_package(CGAL)
if (CGAL_FOUND)
  message("CGAL FOUND.")
  message("CGAL INCLUDE DIR = ${CGAL_INCLUDE_DIRS}")
  message("CGAL LIBRARIES = ${CGAL_LIBRARY}")
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
add_subdirectory(test)
#add_subdirectory(ModelViewer)

add_executable(InstantaneousInterface main.cpp)
add_executable(mesh_op mesh_op.cpp)
#add_executable(mesh_display mesh_display.cpp)

target_include_directories(InstantaneousInterface PUBLIC ${CMAKE_SOURCE_DIR})
target_include_directories(mesh_op PUBLIC ${CMAKE_SOURCE_DIR})
#target_include_directories(mesh_display PUBLIC ${CMAKE_SOURCE_DIR} ${OPENGL_INCLUDE_DIR} ${GLFW3_INCLUDE_DIR})

target_link_libraries(InstantaneousInterface PUBLIC IIxdr IIsrc IItools IIparallel)
target_link_libraries(mesh_op PUBLIC IIxdr IIsrc IItools IIparallel)
#target_link_libraries(mesh_display PUBLIC IIxdr IIsrc IItools IIparallel IIModelViewer ${GLFW3_LIBRARY})


install(TARGETS InstantaneousInterface DESTINATION bin) 
install(TARGETS mesh_op DESTINATION bin)
#install(TARGETS mesh_display DESTINATION bin)
