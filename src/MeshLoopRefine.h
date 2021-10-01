#pragma once
#include "MeshRefineStrategy.h"
#include "Mesh.h"

#include <vector>
#include <string>
#include <cmath>

class MeshLoopRefine : public MeshRefineStrategy
{
    public:
        MeshLoopRefine(MeshRefineStrategyInput& input);

        virtual void refine();

        void calculateEvenPoints();

        Real calculateBeta(int k);
    
    private:
        std::vector<Real3> EvenPoints_;
};