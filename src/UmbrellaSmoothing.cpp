#include "UmbrellaSmoothing.h"

namespace MeshRefineStrategyFactory
{
    registry_<UmbrellaSmoothing> registerUmbrella("umbrella");
}

UmbrellaSmoothing::UmbrellaSmoothing(MeshRefineStrategyInput& input)
:MeshRefineStrategy(input)
{
    input.pack.ReadNumber("iterations", ParameterPack::KeyType::Required, numIterations_);
    input.pack.ReadNumber("lambdadt", ParameterPack::KeyType::Optional, lambdadt_);
}

void UmbrellaSmoothing::refine()
{
    mesh_.findVertexNeighbors();
    mesh_.findBoundaryVertices();
    const auto& vertices = mesh_.getvertices();

    // fill the old vertices with the original vertices
    oldVertices_.insert(oldVertices_.end(), vertices.begin(), vertices.end());

    // resize the new vertices 
    newVertices_.resize(vertices.size());

    for (int i=0;i<numIterations_;i++)
    {
        refineStep();
    }

    auto& vert = mesh_.accessvertices();
    vert.clear();
    vert.insert(vert.end(), newVertices_.begin(), newVertices_.end());
}

void UmbrellaSmoothing::refineStep()
{
    // get the vertices 
    const auto& vertices = mesh_.getvertices();

    // get neighbors indices 
    const auto& neighborIndices = mesh_.getNeighborIndices();

    // clear new vertices 
    newVertices_.clear();
    newVertices_.insert(newVertices_.end(),oldVertices_.begin(), oldVertices_.end());

    #pragma omp parallel for
    for (int i=0;i<vertices.size();i++)
    {
        if (! mesh_.isBoundary(i))
        {
            int numNeighbors = neighborIndices[i].size();

            // the new vertex 
            Real3 newV;
            newV.fill(0);

            Real factor = lambdadt_/numNeighbors;
            Real3 diff;
            diff.fill(0);

            // perform 1/m \sum_{neighbors} xj - xi
            for (int j=0;j<numNeighbors;j++)
            {
                int id = neighborIndices[i][j];
                for (int k=0;k<3;k++)
                {
                    diff[k] += factor * (oldVertices_[id].position_[k] - oldVertices_[i].position_[k]);
                }
            }

            for (int j=0;j<3;j++)
            {
                newV[j] += diff[j] + oldVertices_[i].position_[j];
            }

            newVertices_[i].position_ = newV;
        }
    }

    oldVertices_.clear();
    oldVertices_.insert(oldVertices_.end(), newVertices_.begin(), newVertices_.end());
}