#pragma once
#include "MeshRefineStrategy.h"
#include "Mesh.h"

#include <vector>
#include <string>
#include <cmath>
#include <unordered_map>

class MeshLoopRefine : public MeshRefineStrategy
{
    public:
        MeshLoopRefine(MeshRefineStrategyInput& input);

        virtual void refine();

        void calculateEvenPoints();
        void calculateOddPoints();
        void calculateTriangles();
        void updateMesh();

        Real calculateBeta(int k);
    
    private:
        // the new vertices in the system 
        std::vector<vertex> Points_;

        // the triangles in the system
        std::vector<triangle> triangles_;

        // create a map to map from edges to the new points
        std::unordered_map<edge, vertex> MapEdgesToVertex_;
};