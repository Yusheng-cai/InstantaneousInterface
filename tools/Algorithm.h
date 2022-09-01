#pragma once

#include "CommonTypes.h"

#include <algorithm>
#include <numeric>
#include <vector>
#include <random>

namespace Algorithm
{
    using Real = CommonTypes::Real;
    void Permutation(int max, int numSamples, int numTimes, std::vector<std::vector<int>>& samples);
    void Permutation(int max, int numSamples, std::vector<int>& samples);

    template <std::size_t dim>
    int argmin(std::array<Real,dim>& arr);

    template <std::size_t dim>
    int argmax(std::array<Real,dim>& arr);

    template <typename T>
    int argmin(std::vector<T>& vec);

    template <typename T>
    int argmax(std::vector<T>& vec);

    template <typename T> 
    std::vector<T> arange(T min, T max, T step);
};

template<typename T>
std::vector<T> Algorithm::arange(T start, T stop, T step) {
    std::vector<T> values;
    for (T value = start; value < stop; value += step)
        values.push_back(value);
    return values;
}

template <std::size_t dim>
int Algorithm::argmin(std::array<Real,dim>& arr)
{
    std::array<Real,dim>::iterator it = std::min_element(arr.begin(), arr.end());

    return it - arr.begin();
}

template <std::size_t dim>
int Algorithm::argmax(std::array<Real,dim>& arr)
{
    std::array<Real,dim>::iterator it = std::max_element(arr.begin(), arr.end());

    return it - arr.begin();
}

template <typename T>
int Algorithm::argmax(std::vector<T>& vec)
{
    std::vector<T>::iterator it = std::max_element(vec.begin(), vec.end());

    return it - vec.begin();
}

template <typename T>
int Algorithm::argmin(std::vector<T>& vec)
{
    std::vector<T>::iterator it = std::min_element(vec.begin(), vec.end());

    return it - vec.begin();
}