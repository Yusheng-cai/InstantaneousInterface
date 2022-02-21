#pragma once

#include "MeshRefineStrategy.h"
#include "tools/Assert.h"
#include "tools/CommonTypes.h"
#include "Mesh.h"
#include "Eigen/Core"
#include "Eigen/Sparse"
#include "Eigen/QR"
#include "Eigen/LU"
#include "Eigen/SparseLU"
#include "parallel/OpenMP_buffer.h"
#include "Eigen/SparseCholesky"

#include <vector>
#include <array>
#include <iterator>
#include <string>
#include <cmath>
#include <math.h>

class MeshCurvatureflow : public MeshRefineStrategy
{
    public:
        using triplet = Eigen::Triplet<Real,int>;
        using Sparse_mat = Eigen::SparseMatrix<Real>;
        using Sparse_Chol= Eigen::SimplicialLLT<Sparse_mat>;
        using Sparse_LU = Eigen::SparseLU<Sparse_mat>;

        MeshCurvatureflow(MeshRefineStrategyInput& input);

        virtual void refine(Mesh& mesh) override;

        void refineExplicitStep();
        void refineImplicitStep();

        void getImplicitMatrix();

        // calculate weights of its 1-ring neighbors of vertex i
        std::vector<Real> calculateWeights(int i, std::vector<int>& neighborId, Real3& Lfactor, bool& flag);

    private:
        int numIterations_;
        Real lambdadt_=0.5;

        std::vector<vertex> newVertices_;
        std::vector<Real> TotalArea_;

        std::string solverName_="explicit";

        // eigen triplet is a data structure that is useful for sparse matrix storage (i,j,value)
        std::vector<triplet> triplets_;
        OpenMP::OpenMP_buffer<std::vector<triplet>> triplets_buffer_;

        Eigen::SparseMatrix<Real> L_;

        // initialize the solver
        Eigen::BiCGSTAB<Sparse_mat> solver_;

        // initial volume of the object 
        Real initialVolume_;

        // whether or not we scale the object by volume every step
        bool scale_=false;

        // The curvature we are setting the surface to be 
        Real k0_ = 0;

        // Lfactors 
        std::vector<Real3> Lfactors_;

        // vertex positions --> includes virtual vertices
        std::vector<Real3> vertexPos_;

        // whether or not we are using virtual sites on the boundary
        bool virtualSite_=false;
};