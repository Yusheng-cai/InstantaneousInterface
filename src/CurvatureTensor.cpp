#include "CurvatureTensor.h"

namespace CurvatureRegistry
{
   registry<CurvatureTensor> registerCurvatureTensor("tensor"); 
}

CurvatureTensor::CurvatureTensor(CurvatureInput& input)
:Curvature(input)
{}

void CurvatureTensor::calculate()
{
    mesh_.CalcPerVertexDir();

    const auto& triangles = mesh_.gettriangles();
    const auto& vertices  = mesh_.getvertices();

    curvatureTensor_.resize(triangles.size());
    curvatureVec_.resize(triangles.size());

    for (int i=0;i<triangles.size();i++)
    {
        auto& t = triangles[i].triangleindices_; 
        std::array<Real3,3> edges;

        for (int j=0;j<3;j++)
        {
            edges[0][j] = vertices[t[2]].position_[j] - vertices[t[1]].position_[j];
            edges[1][j] = vertices[t[0]].position_[j] - vertices[t[2]].position_[j];
            edges[2][j] = vertices[t[1]].position_[j] - vertices[t[0]].position_[j];
        }

        Real3 U = edges[0];
        LinAlg3x3::normalize(U);
        Real3 N = LinAlg3x3::CrossProduct(edges[0], edges[1]);
        LinAlg3x3::normalize(N);
        Real3 V = LinAlg3x3::CrossProduct(U, N);
        LinAlg3x3::normalize(V); 

        // initialize A and B matrix to be solved
        Eigen::Matrix3d A;
        Eigen::Vector3d b;
        b.fill(0);
        A.fill({});
 
        for (int j=0;j<3;j++)
        {
            Real ejU = LinAlg3x3::DotProduct(edges[j], U);
            Real ejV = LinAlg3x3::DotProduct(edges[j], V); 
            A(0,0) += ejU*ejU;
            A(0,1) += ejU*ejV;
            A(2,2) += ejV*ejV;

            int id1 = j - 1;
            if ( id1 < 0)
            {
                id1 += 3;
            }

            int id2 = j - 2;
            if ( id2 < 0)
            {
                id2 += 3;
            }
            Real3 diffN;

            for (int k =0;k<3;k++)
            {
                diffN[k] = vertices[t[id1]].normals_[k] - vertices[t[id2]].normals_[k];
            }

            Real diffNU = LinAlg3x3::DotProduct(diffN, U);
            Real diffNV = LinAlg3x3::DotProduct(diffN, V);

            b[0] += diffNU*ejU;
            b[1] += diffNU*ejV + diffNV*ejU;
            b[2] += diffNV*ejV;
        }

        A(1,1) = A(0,0) + A(2,2);
        A(1,2) = A(0,1);
        A(2,1) = A(0,1);
        A(1,0) = A(0,1);

        Eigen::Vector3d soln = A.bdcSvd(Eigen::ComputeThinU|Eigen::ComputeThinV).solve(b);
        Real3 ans;

        for (int j=0;j<3;j++)
        {
            ans[j] = soln[j];
        }

        curvatureTensor_[i] = ans;
 
        Eigen::Matrix2d mat;
        Eigen::EigenSolver<Eigen::Matrix2d> eigensolver;

        mat(0,0) = ans[0];
        mat(0,1) = ans[1];
        mat(1,0) = ans[1];
        mat(1,1) = ans[2];
        eigensolver.compute(mat);

        Eigen::Vector2d eigenvalues = eigensolver.eigenvalues().real();

        curvatureVec_[i][0] = eigenvalues[0];
        curvatureVec_[i][1] = eigenvalues[1];
    }
}

void CurvatureTensor::printOutput()
{
    if (ofs_.is_open())
    {
        for (int i=0;i<curvatureVec_.size();i++)
        {
            for (int j=0;j<2;j++)
            {
                ofs_ << curvatureVec_[i][j] << " ";
            }
            ofs_ << "\n";
        }
        ofs_.close();
    }
}