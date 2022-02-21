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
        using index3 = CommonTypes::index3;
        using index2 = CommonTypes::index2;
        MeshLoopRefine(MeshRefineStrategyInput& input);

        virtual void refine(Mesh& mesh);

        void calculateEvenPoints();
        void calculateOddPoints();
        void calculateTriangles();
        void updateMesh();

        Real calculateBeta(int k);
    
    private:
        // the new vertices in the system 
        std::vector<vertex> Points_;

        // nullptr for mesh obj
        Mesh* mesh_ = nullptr;

        // the triangles in the system
        std::vector<triangle> triangles_;
        std::vector<index3> triangle_indices_;

        // create a map to map from edges to the new points
        std::map<index2, vertex> MapEdgesToVertex_;
};