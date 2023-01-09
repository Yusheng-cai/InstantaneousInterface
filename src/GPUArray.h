#pragma once

#include "tools/Assert.h"

#include <iostream>
#include <vector>
#include <cuda_runtime.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>

template <typename T>
class GPUArray2d{
    public:
        using INT2 = std::array<int,2>;
        GPUArray2d(INT2 sizes);
        GPUArray2d(int size1, int size2);
        GPUArray2d(INT2 sizes, T num);
        GPUArray2d(int size1, int size2, T num);

        T& operator()(int index1, int index2);
        T& operator()(INT2 indices);
        const T& operator()(int index1, int index2) const;
        const T& operator()(INT2 indices) const;

        void memcpyHostToDevice();
        void memcpyDeviceToHost();

        T* host_data();
        T* device_data();

        INT2 getSize() {return {{size1_, size2_}};}
        int getTotalSize() {return total_size_;}

        int CalculateIndex(int index1, int index2);

    private:
        int size1_, size2_;
        int total_size_;
        thrust::device_vector<T> device_vector_;
        thrust::host_vector<T> host_vector_;
};


                        /// 2d  array ///

template <typename T>
GPUArray2d<T>::GPUArray2d(INT2 sizes)
:size1_(sizes[0]), size2_(sizes[1])
{
    // calculate the total size of the array
    total_size_ = size1_ * size2_;
    host_vector_.resize(total_size_);
}

template <typename T>
GPUArray2d<T>::GPUArray2d(int size1, int size2)
:size1_(size1), size2_(size2)
{
    // calculate the total size of the array
    total_size_ = size1_ * size2_;
    host_vector_.resize(total_size_);
}


template <typename T>
GPUArray2d<T>::GPUArray2d(INT2 sizes, T num)
:GPUArray2d(sizes)
{
    for (int i=0;i<total_size_;i++){
        host_vector_[i] = num;
    }
}

template <typename T>
GPUArray2d<T>::GPUArray2d(int size1, int size2, T num)
:GPUArray2d(size1,size2)
{
    for (int i=0;i<total_size_;i++){
        host_vector_[i] = num;
    }
}

template <typename T>
void GPUArray2d<T>::memcpyHostToDevice(){
    device_vector_ = host_vector_;
}

template <typename T>
void GPUArray2d<T>::memcpyDeviceToHost(){
    thrust::copy(device_vector_.begin(), device_vector_.end(), host_vector_.begin());
}

template <typename T>
T& GPUArray2d<T>::operator()(int index1, int index2){
    int ind = CalculateIndex(index1, index2);

    return host_vector_[ind];
}

template <typename T>
T& GPUArray2d<T>::operator()(INT2 indices){
    return GPUArray2d<T>::operator()(indices[0], indices[1]);
}

template <typename T>
const T& GPUArray2d<T>::operator()(int index1, int index2) const {
    int ind = CalculateIndex(index1, index2);

    return host_vector_[ind];
}

template <typename T>
const T& GPUArray2d<T>::operator()(INT2 indices) const{
    return GPUArray2d<T>::operator()(indices[0], indices[1]);
}


template <typename T>
int GPUArray2d<T>::CalculateIndex(int index1, int index2){
    int index1d = index1 * size2_ + index2;
    ASSERT((index1d < total_size_), "Index " << index1 << " " << index2 <<  " out of range");

    return index1d;
}

template <typename T>
T* GPUArray2d<T>::host_data(){
    return thrust::raw_pointer_cast(host_vector_.data());
}

template <typename T>
T* GPUArray2d<T>::device_data(){
    return thrust::raw_pointer_cast(device_vector_.data());
}