#include "MeshLoopRefine.h"

MeshLoopRefine::MeshLoopRefine(MeshRefineStrategyInput& input)
:MeshRefineStrategy(input)
{
}

void MeshLoopRefine::refine()
{
}

void MeshLoopRefine::calculateEvenPoints()
{
    // first we find the neighbors of the vertex
    mesh_.findVertexNeighbors();

    const auto& neighbors = mesh_.getNeighborIndices();
    const auto& vertices  = mesh_.getvertices();

    // resize the evenpoints
    int numPoints = neighbors.size();
    EvenPoints_.resize(numPoints,{});

    for (int i=0;i<numPoints;i++)
    {
        int lenNeighbors = neighbors[i].size();
        Real beta = this->calculateBeta(lenNeighbors);

        Real3 posI = vertices[i].position_;
        Real3 sum_ = {};
        for (int j=0;j<lenNeighbors;j++)
        {
            Real3 posN = vertices[neighbors[i][j]].position_;
            for (int k=0;k<3;k++)
            {
                sum_[k] += posN[k];
            }
        }

        for (int j=0;j<3;j++)
        {
            sum_[j] = (1 - lenNeighbors*beta)*posI[j] + beta * sum_[j];
        }
        EvenPoints_[i] = sum_;
    }
}

MeshLoopRefine::Real MeshLoopRefine::calculateBeta(int n)
{
    if (n == 3)
    {
        return 3.0/16.0;
    }
    else
    {
        Real factor = 1.0/(Real)n;
        Real l = 3.0/8.0 + 0.25*std::cos(6.28318531/(Real)n);
        Real l2 = l*l;

        Real beta = factor * (5.0/8.0 - 0.25*l2);

        return beta;
    }
}
