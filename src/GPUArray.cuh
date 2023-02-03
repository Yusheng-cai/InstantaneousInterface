#pragma once

#include "tools/Assert.h"

#include <iostream>
#include <vector>
#include <cuda_runtime.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/sort.h>
#include <thrust/execution_policy.h>

template <typename T>
class GPUArray1d
{
    public:
        GPUArray1d() {size_=0;}
        GPUArray1d(int size);

        T& operator()(int index1);
        const T& operator()(int index1) const;

        void memcpyHostToDevice();
        void memcpyDeviceToHost();

        void resize(int index1);

        T* host_data();
        T* device_data();
        thrust::device_vector<T>& device_vector() {return device_vector_;}
        thrust::host_vector<T>& host_vector() {return host_vector_;}

        int getTotalSize() {return size_;}

    private:
        int size_;
        thrust::device_vector<T> device_vector_;
        thrust::host_vector<T> host_vector_;
};

template <typename T>
GPUArray1d<T>::GPUArray1d(int size)
:size_(size){
    host_vector_.resize(size_);
}

template <typename T>
T& GPUArray1d<T>::operator() (int index){
    return host_vector_[index];
}

template <typename T>
const T& GPUArray1d<T>::operator() (int index) const{
    return host_vector_[index];
}

template <typename T>
void GPUArray1d<T>::resize(int size){
    host_vector_.resize(size);
    size_ = size;
}

template <typename T>
T* GPUArray1d<T>::host_data(){
    return thrust::raw_pointer_cast(host_vector_.data());
}

template <typename T>
T* GPUArray1d<T>::device_data(){
    return thrust::raw_pointer_cast(device_vector_.data());
}

template <typename T>
void GPUArray1d<T>::memcpyHostToDevice(){
    device_vector_ = host_vector_;
}

template <typename T>
void GPUArray1d<T>::memcpyDeviceToHost(){
    thrust::copy(device_vector_.begin(), device_vector_.end(), host_vector_.begin());
}


template <typename T>
class GPUArray2d
{
    public:
        using INT2 = std::array<int,2>;
        GPUArray2d() {total_size_=0;}
        GPUArray2d(INT2 sizes);
        GPUArray2d(int size1, int size2);
        GPUArray2d(INT2 sizes, T num);
        GPUArray2d(int size1, int size2, T num);

        T& operator()(int index1, int index2);
        T& operator()(INT2 indices);
        const T& operator()(int index1, int index2) const;
        const T& operator()(INT2 indices) const;

        T& operator[](int index) {return host_vector_[index];}
        const T& operator[](int index) const {return host_vector_[index];}

        void memcpyHostToDevice();
        void memcpyDeviceToHost();

        void resize(int index1, int index2);
        void resize(INT2 indices);

        T* host_data();
        T* device_data();
        thrust::device_vector<T>& device_vector() {return device_vector_;}
        thrust::host_vector<T>& host_vector() {return host_vector_;}

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
void GPUArray2d<T>::resize(int size1, int size2){
    size1_ = size1; size2_=size2;
    total_size_=size1_ * size2_;
    host_vector_.resize(total_size_);
}

template <typename T>
void GPUArray2d<T>::resize(INT2 indices){
    resize(indices[0], indices[1]);
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

/////// 3d 
template <typename T>
class GPUArray3d
{
    public:
        using INT3 = std::array<int,3>;
        GPUArray3d() {total_size_=0;}
        GPUArray3d(INT3 sizes);
        GPUArray3d(int size1, int size2, int size3);
        GPUArray3d(INT3 sizes, T num);
        GPUArray3d(int size1, int size2, int size3, T num);

        T& operator()(int index1, int index2, int index3);
        T& operator()(INT3 indices);
        T& operator[](int ind) {return host_vector_[ind];}
        const T& operator[](int ind) const {return host_vector_[ind];}

        const T& operator()(int index1, int index2, int index3) const;
        const T& operator()(INT3 indices) const;

        void resize(int index1, int index2, int index3);
        void resize(int index1, int index2, int index3, T num);
        void resize(INT3 indices);

        void memcpyHostToDevice();
        void memcpyDeviceToHost();

        T* host_data();
        T* device_data();
        thrust::device_vector<T>& device_vector() {return device_vector_;}
        thrust::host_vector<T>& host_vector() {return host_vector_;}

        INT3 getSize() {return {{size1_, size2_, size3_}};}
        int getTotalSize() {return total_size_;}

        int CalculateIndex(int index1, int index2, int index3);

    private:
        int size1_, size2_, size3_;
        int total_size_;
        thrust::device_vector<T> device_vector_;
        thrust::host_vector<T> host_vector_;
};

template <typename T>
GPUArray3d<T>::GPUArray3d(INT3 sizes)
:size1_(sizes[0]), size2_(sizes[1]), size3_(sizes[2])
{
    // calculate the total size of the array
    total_size_ = size1_ * size2_ * size3_;
    host_vector_.resize(total_size_);
}


template <typename T>
GPUArray3d<T>::GPUArray3d(int size1, int size2, int size3)
:size1_(size1), size2_(size2), size3_(size3)
{
    // calculate the total size of the array
    total_size_ = size1_ * size2_ * size3_;
    host_vector_.resize(total_size_);
}


template <typename T>
GPUArray3d<T>::GPUArray3d(INT3 sizes, T num)
:GPUArray3d(sizes)
{
    for (int i=0;i<total_size_;i++){
        host_vector_[i] = num;
    }
}

template <typename T>
GPUArray3d<T>::GPUArray3d(int size1, int size2, int size3, T num)
:GPUArray3d(size1,size2,size3)
{
    for (int i=0;i<total_size_;i++){
        host_vector_[i] = num;
    }
}

template <typename T>
void GPUArray3d<T>::resize(int size1, int size2, int size3){
    size1_=size1; size2_=size2; size3_=size3;
    total_size_=size1_*size2_*size3_;
    host_vector_.resize(total_size_);
}

template <typename T>
void GPUArray3d<T>::resize(int size1, int size2, int size3, T num){
    size1_=size1; size2_=size2; size3_=size3;
    total_size_=size1_*size2_*size3_;
    host_vector_.resize(total_size_, num);
}

template <typename T>
void GPUArray3d<T>::resize(INT3 indices){
    resize(indices[0], indices[1], indices[2]);
}

template <typename T>
void GPUArray3d<T>::memcpyHostToDevice(){
    device_vector_ = host_vector_;
}

template <typename T>
void GPUArray3d<T>::memcpyDeviceToHost(){
    thrust::copy(device_vector_.begin(), device_vector_.end(), host_vector_.begin());
}

template <typename T>
T& GPUArray3d<T>::operator()(int index1, int index2, int index3){
    int ind = CalculateIndex(index1, index2, index3);

    return host_vector_[ind];
}

template <typename T>
T& GPUArray3d<T>::operator()(INT3 indices){
    return GPUArray3d<T>::operator()(indices[0], indices[1], indices[2]);
}

template <typename T>
const T& GPUArray3d<T>::operator()(int index1, int index2, int index3) const {
    int ind = CalculateIndex(index1, index2, index3);

    return host_vector_[ind];
}

template <typename T>
const T& GPUArray3d<T>::operator()(INT3 indices) const{
    return GPUArray3d<T>::operator()(indices[0], indices[1], indices[2]);
}


template <typename T>
int GPUArray3d<T>::CalculateIndex(int index1, int index2, int index3){
    int index1d = index3 * size1_ * size2_ + index2 * size1_ + index1;
    ASSERT((index1d < total_size_), "Index " << index1 << " " << index2  << " " << index3 <<  " out of range");

    return index1d;
}

template <typename T>
T* GPUArray3d<T>::host_data(){
    return thrust::raw_pointer_cast(host_vector_.data());
}

template <typename T>
T* GPUArray3d<T>::device_data(){
    return thrust::raw_pointer_cast(device_vector_.data());
}
