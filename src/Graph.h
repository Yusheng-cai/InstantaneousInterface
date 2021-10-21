#pragma once
#include "tools/Assert.h"
#include "tools/CommonTypes.h"

#include <cmath>
#include <iostream>
#include <vector>
#include <string>
#include <algorithm>

//#include "DistanceMatrix.h"

namespace Graph
{
    using index3 = CommonTypes::index3;
    // This assumes that the vector is shaped accordingly
    // template<typename T>
    // void FloryWarshallAlgo(const PairDistanceMatrix<T>& distance, PairDistanceMatrix<T>& graph);

    // Need the input of nearby Indices of the vertex connected to it 
    void getNearbyIndicesNVertexAway(const std::vector<std::vector<int>>& NearbyIndices, int N, std::vector<std::vector<int>>& NerabyIndicesNVertex);

    void getNNearbyIndices(const std::vector<std::vector<int>>& NearbyIndices, int N, std::vector<std::vector<int>>& NearbyIndicesN);
};

// template<typename T>
// void Graph::FloryWarshallAlgo(const PairDistanceMatrix<T>& distance, PairDistanceMatrix<T>& graph)
// {   
//     int N = distance.getSize();
//     graph = PairDistanceMatrix<T>(distance);

//     double INF = PairDistanceMatrix<T>::INF;

//     // iterate over the vertices
//     for (int i=0;i<N;i++)
//     {
//         // iterate over the start
//         for (int j=0;j<N;j++)
//         {
//             // iterate over the end
//             for (int k=0;k<N;k++)
//             {
//                 if ((graph(j,k) > (graph(j,i) + graph(i,k))) && graph(j,i) != INF && graph(i,k) != INF)
//                 {
//                     graph(i,j) = graph(i,k) + graph(k,j);
//                 }
//             }
//         }
//     }
// }