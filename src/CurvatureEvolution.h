#pragma once 

#include "MeshRefineStrategy.h"
#include "Curvature.h"
#include "tools/CommonTypes.h"
#include "xdr/XtcFile.h"

#include <vector>
#include <string>
#include <array>
#include <memory>
#include <chrono>
#include <algorithm>
#include <numeric>

// We are trying to solve the general equation of 
// dxi/dt = (k0 - kmean) * ni
// However, our kmean could be what we calculate 

class CurvatureEvolution : public MeshRefineStrategy
{
    public:
        using curveptr = std::unique_ptr<Curvature>;
        using xtcptr   = std::unique_ptr<XtcFile>;
        using Matrix   = CommonTypes::Matrix;
        using Real2    = CommonTypes::Real2;

        CurvatureEvolution(MeshRefineStrategyInput& input);

        virtual void refine(Mesh& mesh) override;

        void findVertices();

    private:
        std::string curvaturetype_;
        curveptr curvatureCalc_;

        // The indices for the outside vertices to not be kept constant 
        std::vector<int> VertexIndices_;

        Real StepSize_;
        Real meanCurvature_;

        Real tol_ = 1e-4;

        // keep track of number of iteration
        int iteration_=0;
        int maxStep = 1e6;

        // whether or not we want to print the error 
        Real err_;
        int skip_=1;
};