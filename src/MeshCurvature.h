#pragma once 

#include "MeshRefineStrategy.h"
#include "Curvature.h"
#include "tools/CommonTypes.h"
#include "xdr/XtcFile.h"

#include <vector>
#include <string>
#include <array>
#include <memory>

// We are trying to solve the general equation of 
// dxi/dt = (k0 - kmean) * ni
// However, our kmean could be what we calculate 

class MeshCurvature : public MeshRefineStrategy
{
    public:
        using curveptr = std::unique_ptr<Curvature>;
        using xtcptr   = std::unique_ptr<XtcFile>;
        using Matrix   = CommonTypes::Matrix;

        MeshCurvature(MeshRefineStrategyInput& input);

        virtual void refine() override;
    private:
        std::string curvaturetype_;
        curveptr curvatureCalc_;

        Real StepSize_;
        Real meanCurvature_;

        Real tol_ = 1e-4;

        std::string xtcName="test.xtc";
        bool xtcWrite_=false;
        xtcptr xtcFile_;
        int xtcstep_=1;
        int skip_=100;
        int maxStep = 1e6;

        Matrix box={{{10,0,0},{0,10,0},{0,0,10}}};
};