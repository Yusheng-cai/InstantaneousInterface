#include "Graph.h"

void Graph::getNearbyIndicesNVertexAway(const std::vector<std::vector<int>>& NearbyIndices, int N, std::vector<std::vector<int>>& NearbyIndicesNVertex)
{
    int NumberOfVertices = NearbyIndices.size();

    NearbyIndicesNVertex.clear();
    NearbyIndicesNVertex.resize(NumberOfVertices);

    #pragma omp parallel for
    for(int i=0;i<NumberOfVertices;i++)
    {
        int initial_count = 0;
        int final_count = 0;
        int nearbySize = NearbyIndices[i].size();
        const std::vector<int>& currentNearbyIndices = NearbyIndices[i];
        std::vector<int> currentNearbyIndicesNVertex;
        currentNearbyIndicesNVertex.resize(0);

        for (int j=0;j<N;j++)
        {
            if (currentNearbyIndicesNVertex.size() == 0)
            {
                currentNearbyIndicesNVertex.insert(currentNearbyIndicesNVertex.end(),\
                currentNearbyIndices.begin(), currentNearbyIndices.end());

                final_count = currentNearbyIndicesNVertex.size();
            }            
            else
            {
                std::vector<int> temp;
                temp.resize(0);
                for (int k=initial_count;k<final_count;k++)
                {
                    int id = currentNearbyIndicesNVertex[k];

                    for (int m=0;m<NearbyIndices[id].size();m++)
                    {
                        if (NearbyIndices[id][m] != i)
                        {
                            temp.push_back(NearbyIndices[id][m]);
                        }
                    }
                }

                currentNearbyIndicesNVertex.insert(currentNearbyIndicesNVertex.end(), temp.begin(), temp.end());
                initial_count = final_count;
                final_count = currentNearbyIndicesNVertex.size();
            }
        }

        std::sort(currentNearbyIndicesNVertex.begin(), currentNearbyIndicesNVertex.end());
        // The non-unique terms will be undefined ***
        auto it = std::unique(currentNearbyIndicesNVertex.begin(), currentNearbyIndicesNVertex.end());

        // removes the non-unique terms
        currentNearbyIndicesNVertex.resize(std::distance(currentNearbyIndicesNVertex.begin(), it));

        NearbyIndicesNVertex[i] = currentNearbyIndicesNVertex;
    }
}