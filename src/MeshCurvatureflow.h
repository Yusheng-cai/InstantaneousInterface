#pragma once

#include "MeshRefineStrategy.h"
#include "tools/Assert.h"
#include "tools/CommonTypes.h"
#include "Mesh.h"


#include <vector>
#include <array>
#include <iterator>
#include <string>

class MeshCurvatureflow : public MeshRefineStrategy
{
    public:
        MeshCurvatureflow(MeshRefineStrategyInput& input);

        virtual void refine() override;

        void refineStep();


    private:
        int numIterations_;
        Real lambdadt_=0.5;

        std::vector<vertex> newVertices_;
        std::vector<Real> TotalArea_;
};
