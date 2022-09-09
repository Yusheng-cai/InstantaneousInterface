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

    template <typename T>
    std::vector<T> linspace(T min, T max, int num);

    template <typename T>
    T max(std::vector<T>& vec);

    template <typename T>
    T min(std::vector<T>& vec);

    template <typename T>
    bool contain(std::vector<T>& vec, T num);

    template <typename T>
    void unique(std::vector<T>& vec); 
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
    typename std::array<Real,dim>::iterator it = std::min_element(arr.begin(), arr.end());

    return it - arr.begin();
}

template <std::size_t dim>
int Algorithm::argmax(std::array<Real,dim>& arr)
{
    typename std::array<Real,dim>::iterator it = std::max_element(arr.begin(), arr.end());

    return it - arr.begin();
}

template <typename T>
int Algorithm::argmax(std::vector<T>& vec)
{
    typename std::vector<T>::iterator it = std::max_element(vec.begin(), vec.end());

    return it - vec.begin();
}

template <typename T>
int Algorithm::argmin(std::vector<T>& vec)
{
    typename std::vector<T>::iterator it = std::min_element(vec.begin(), vec.end());

    return it - vec.begin();
}

template <typename T>
std::vector<T> Algorithm::linspace(T min, T max, int num)
{
    T step = (max - min)/num;
    std::vector<T> ret;
    for (int i=0;i<num;i++)
    {
        ret.push_back(min + i*step);
    }

    return ret;
}

template <typename T>
T Algorithm::max(std::vector<T>& vec)
{
    typename std::vector<T>::iterator max_element = std::max_element(vec.begin(), vec.end());

    return *max_element;
}

template <typename T>
T Algorithm::min(std::vector<T>& vec)
{
    typename std::vector<T>::iterator min_element = std::min_element(vec.begin(), vec.end());

    return *min_element;
}

template <typename T>
bool Algorithm::contain(std::vector<T>& vec, T num)
{
    return (std::find(vec.begin(), vec.end(), num) != vec.end());
}

template <typename T>
void Algorithm::unique(std::vector<T>& vec)
{
    std::sort(vec.begin(), vec.end());
    typename std::vector<T>::iterator it = std::unique(vec.begin(), vec.end());
    vec.resize(std::distance(vec.begin(), it));
}