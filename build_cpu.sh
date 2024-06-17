#!/bin/bash

export CC=gcc
export CXX=g++

module load cuda/11.7

# specify the build directory
build_type=RELEASE
build_dir=$PWD/${build_type}_cpu/
install_dir=${HOME}/programs/InstantaneousInterface/cpu
libigl_dir=${HOME}/programs/libigl/

# remove build_dir if it already exists
[[ -d ${build_dir} ]] && rm -rf ${build_dir}
mkdir -p $build_dir

cd $build_dir
cmake .. -DCMAKE_BUILD_TYPE=${build_type} \
	 -DCMAKE_INSTALL_PREFIX=${install_dir} \
	 #-DLIBIGL_DIR=${libigl_dir} \
	 #-DCMAKE_CUDA_COMPILER="/usr/local/cuda-11.7/bin/nvcc" \
	 -DENABLE_CUDA=OFF \
	 -DENABLE_IGL=OFF
	 #-DGLFW3_ROOT="/home/yusheng/programs/glfw-3.3.8/"

# make with 8 threads
make  -j 24

make test
make install
