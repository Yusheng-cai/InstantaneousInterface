#pragma once 

#include "MeshRefineStrategy.h"
#include "Curvature.h"
#include "tools/CommonTypes.h"
#include "xdr/XtcFile.h"
#include "MeshCurvatureflow.h"
#include "ShortEdgeRemoval.h"
#include "ObtuseTriangleRemoval.h"

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

class InterfacialFE_minimization : public MeshRefineStrategy{
    public:
        InterfacialFE_minimization(MeshRefineStrategyInput& input);

        void update_Mesh();

        virtual void refine(Mesh& mesh) override;

        const std::vector<Real>& getFE() {return FE_;}

        // function to obtain the gamma
        Real CalculateGamma(Real temperature);
        Real CalculateMu(Real temperature);

        Real getL() {return L_;}
        void setL(Real L) {L_=L;}

        Real getrho() {return rho_;}
        Real getmu() {return mu_;}
        Real getgamma() {return gamma_;}


    private:
        // triangle areas 
        std::vector<Real> TriangleAreas_;
        std::vector<std::vector<int>> MapVertexToFace_;
        std::map<INT2, std::vector<int>> MapEdgeToOpposingVerts_;

        int numVerts_, maxStep_=1e5;
        Real tol_=1e-5;
        Real stepsize_=1e-1, k0_;

        std::vector<vertex> newVertices_;
        std::vector<Real> TotalArea_;

        std::vector<Real3> dAdpi_, dVdpi_;

        Real temperature_=298.0f;
        Real L_=0.0f;
        Real gamma_; // kJ/nm2
        Real mu_;
        Real rho_=5.0333e-23; // mol/nm3

        int print_every=100, optimize_every=1e10;

        std::vector<Real> FE_;

        bool MaxStepCriteria=true;
};