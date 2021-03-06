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

        void refineImplicitStep();

        Eigen::SparseMatrix<Real> CalculateImplicitMatrix();

        // calculate weights of its 1-ring neighbors of vertex i
        bool CalculateCotangentWeights(int i, Real3& Lfactor, std::vector<Real>& cotFactors);

    private:
        int numIterations_;
        Real lambdadt_=0.5;

        std::vector<vertex> newVertices_;
        std::vector<Real> TotalArea_;

        // eigen triplet is a data structure that is useful for sparse matrix storage (i,j,value)
        std::vector<triplet> triplets_;

        // The sparse matrix L 
        Eigen::SparseMatrix<Real> L_;

        // initial volume of the object 
        Real initialVolume_;

        // whether or not we scale the object by volume every step
        bool scale_=false;

        // Lfactors 
        std::vector<Real3> Lfactors_;

        // triangle areas 
        std::vector<Real> TriangleAreas_;
        std::vector<std::vector<int>> MapVertexToFace_;
        std::map<INT2, std::vector<int>> MapEdgeToOpposingVerts_;

        // number of vertices and faces
        int numVerts_, numFaces_;

        // factor at which we detect INF
        Real epsilon_=1e-8;

        // rhs of the equation as well as the solved results 
        std::vector<Eigen::VectorXf> rhs_;
        std::vector<Eigen::VectorXf> xyz_;
};