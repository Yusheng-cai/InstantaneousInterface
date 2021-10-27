#pragma once

#include "MeshRefineStrategy.h"
#include "tools/Assert.h"
#include "tools/CommonTypes.h"
#include "Mesh.h"

#include <vector>
#include <array>
#include <string>
#include <iostream>

// umbrella smoothing 
/*
    umbrella smoothing attempts to solve the following differential equation 

    dX/dt = lambda L(X)

    where L is the umbrella operator, namely L(xi) = 1/m \sum_{j=1}^{N1(i)} xj - xi where N1(i) signifies the 
    1-ring neighbors of vertex i

    The iterative scheme is as follow Xn+1 = (I + \lambda dt * L)Xn
*/

class UmbrellaSmoothing : public MeshRefineStrategy
{
    public:
        UmbrellaSmoothing(MeshRefineStrategyInput& input);

        virtual void refine() override;

        void refineStep();
    
    private:
        int numIterations_;
        
        Real lambdadt_=0.5;

        std::vector<vertex> newVertices_;
        std::vector<vertex> oldVertices_;
};