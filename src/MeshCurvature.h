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

// We are trying to solve the general equation of 
// dxi/dt = (k0 - kmean) * ni
// However, our kmean could be what we calculate 

class MeshCurvature : public MeshRefineStrategy
{
    public:
        using curveptr = std::unique_ptr<Curvature>;
        using xtcptr   = std::unique_ptr<XtcFile>;
        using Matrix   = CommonTypes::Matrix;
        using Real2    = CommonTypes::Real2;

        MeshCurvature(MeshRefineStrategyInput& input);

        virtual void refine() override;

        void findVertices();

        void initializeErrorFile();

        void printToErrorFile();

    private:
        std::string curvaturetype_;
        curveptr curvatureCalc_;

        // The indices for the vertices to be kept constant 
        std::vector<int> VertexIndices_;

        // The indices for the outside vertices to not be kept constant 
        std::vector<int> OutsideIndices_;

        Real StepSize_;
        Real meanCurvature_;

        Real tol_ = 1e-4;

        // keep track of number of iteration
        int iteration_=0;

        std::string xtcName="test.xtc";
        bool xtcWrite_=false;
        xtcptr xtcFile_;
        int xtcstep_=1;
        int xtcskip_=100;
        int maxStep = 1e6;
        bool ignoreOppositeSign_=false;
        bool fixBoundary_=true;

        // whether or not we want to print the error 
        std::string errorFile_;
        std::vector<Real> error_;
        Real err_;
        std::ofstream errorStream_;
        int skip_=1;

        Matrix box={{{10,0,0},{0,2,0},{0,0,2}}};
};