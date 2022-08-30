#pragma once

#include <algorithm>
#include <numeric>
#include <vector>
#include <random>

namespace Algorithm
{
    void Permutation(int max, int numSamples, int numTimes, std::vector<std::vector<int>>& samples);
    void Permutation(int max, int numSamples, std::vector<int>& samples);

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