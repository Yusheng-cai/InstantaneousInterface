#pragma once

#include <algorithm>
#include <numeric>
#include <vector>
#include <random>

namespace Algorithm
{
    void Permutation(int max, int numSamples, int numTimes, std::vector<std::vector<int>>& samples);
    void Permulation(int max, int numSamples, std::vector<int>& samples);
};