#pragma once
#include "Graph.h"
#include "Curvature.h"
#include "tools/Assert.h"
#include "tools/CommonTypes.h"
#include "Graph.h"
#include "LinAlgTools.h"

#include <vector>
#include <array>

class CurvatureCurveFit : public Curvature
{
    public:
        using Matrix = CommonTypes::Matrix;
        using Real2  = CommonTypes::Real2;
        CurvatureCurveFit(CurvatureInput& input);

        virtual void calculate() override;
        virtual void printCurvature(std::string name) override;
        void printPrincipalDir(std::string name);
    
    private:
        int NumParameters_ = 3;
        int NumNeighbors_ = 1;

        std::vector<Real2> CurvaturePerVertex_;
        std::vector<Real3> principalDir1_;
        std::vector<Real3> principalDir2_;
};