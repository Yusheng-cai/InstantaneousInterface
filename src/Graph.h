#pragma once
#include "tools/Assert.h"
#include "tools/CommonTypes.h"

#include <cmath>
#include <iostream>
#include <vector>
#include <string>
#include <algorithm>

namespace Graph
{
    using index3 = CommonTypes::index3;
    using Real   = CommonTypes::Real;
    using Real3  = CommonTypes::Real3;

    // Need the input of nearby Indices of the vertex connected to it 
    void getNearbyIndicesNVertexAway(const std::vector<std::vector<int>>& NearbyIndices, int N, std::vector<std::vector<int>>& NerabyIndicesNVertex);

    void getNNearbyIndices(const std::vector<std::vector<int>>& NearbyIndices, int N, std::vector<std::vector<int>>& NearbyIndicesN);
};