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


    private:
        // triangle areas 
        std::vector<Real> TriangleAreas_;
        std::vector<std::vector<int>> MapVertexToFace_;
        std::map<INT2, std::vector<int>> MapEdgeToOpposingVerts_;

        int numVerts_, maxStep_=1e5;
        Real stepsize_=1e-1, k0_;
        

        std::vector<vertex> newVertices_;
        std::vector<Real> TotalArea_;

        std::vector<Real3> dAdpi_, dVdpi_;
};