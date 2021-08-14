#include "Graph.h"

void Graph::getNearbyIndicesNVertexAway(const std::vector<std::vector<int>>& NearbyIndices, int N, std::vector<std::vector<int>>& NearbyIndicesNVertex)
{
    int NumberOfVertices = NearbyIndices.size();

    NearbyIndicesNVertex.clear();
    NearbyIndicesNVertex.resize(NumberOfVertices);

    for(int i=0;i<NumberOfVertices;i++)
    {
        int initial_count = 0;
        int final_count = 0;
        int nearbySize = NearbyIndices[i].size();
        const std::vector<int>& currentNearbyIndices = NearbyIndices[i];
        std::vector<int>& currentNearbyIndicesNVertex = NearbyIndicesNVertex[i];
        currentNearbyIndicesNVertex.resize(0);

        for (int j=0;j<nearbySize;j++)
        {
            if (currentNearbyIndices.size() == 0)
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
                    temp.insert(temp.end(), NearbyIndices[id].begin(), NearbyIndices[id].end());
                }

                currentNearbyIndicesNVertex.insert(currentNearbyIndicesNVertex.end(), temp.begin(), temp.end());
                initial_count = final_count;
                final_count = currentNearbyIndicesNVertex.size();
            }
        }
        // The non-unique terms will be undefined ***
        auto it = std::unique(currentNearbyIndicesNVertex.begin(), currentNearbyIndicesNVertex.end());

        // removes the non-unique terms
        currentNearbyIndicesNVertex.resize(std::distance(currentNearbyIndicesNVertex.begin(), it));
    }
}