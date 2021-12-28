#include "CurvatureCurveFit.h"

namespace CurvatureRegistry
{
    registry<CurvatureCurveFit> registerCurvefit("curvefit");
}

CurvatureCurveFit::CurvatureCurveFit(CurvatureInput& input)
:Curvature(input)
{
    input.pack.ReadNumber("neighbors", ParameterPack::KeyType::Optional, NumNeighbors_);
}

void CurvatureCurveFit::calculate(Mesh& mesh)
{
    initialize(mesh);

    mesh.findVertexNeighbors();

    const auto& VertexNeighbors_ = mesh.getNeighborIndices();
    const auto& vertices = mesh.getvertices();

    std::vector<std::vector<int>> NeighborIndicesNVertex_;
    Graph::getNearbyIndicesNVertexAway(VertexNeighbors_, NumNeighbors_,NeighborIndicesNVertex_);

    // the reference direction of the normal vector is the z vector
    Real3 referenceDir = {{0,0,1}};

    for (int i=0;i<vertices.size();i++)
    {
        auto& v = vertices[i];
        auto& neighbors = NeighborIndicesNVertex_[i];

        ASSERT((neighbors.size()>=3), "The number of points to be fit must be larger or equal to 3.");

        // How to rotate normal vector to be z vector
        Matrix rotationMat = LinAlg3x3::GetRotationMatrix(v.normals_, referenceDir); 
        Matrix revrotMat   = LinAlg3x3::GetRotationMatrix(referenceDir, v.normals_);
        Real3 vector = LinAlg3x3::MatrixDotVector(rotationMat, v.normals_);

        Eigen::Matrix3d mat = Eigen::Matrix3d::Zero();
        Eigen::Vector3d b = Eigen::Vector3d::Zero();

        for (int j=0;j<neighbors.size();j++)
        {
            int neighborId = neighbors[j];
            Real3 neighborPos = vertices[neighborId].position_;
            Real3 diff;
            diff.fill(0);

            for (int k=0;k<3;k++)
            {
                diff[k] = neighborPos[k] - v.position_[k];
            }

            // Rotate neighbor position to the frame of reference of interest
            Real3 neighborRotatedPos = LinAlg3x3::MatrixDotVector(rotationMat, diff);

            mat(0,0) += 0.25*std::pow(neighborRotatedPos[0],4.0);
            mat(0,1) += 0.5*std::pow(neighborRotatedPos[0],3.0)*neighborRotatedPos[1];
            mat(0,2) += 0.25*std::pow(neighborRotatedPos[0],2.0)*std::pow(neighborRotatedPos[1],2.0);
            mat(1,0) += 0.5*std::pow(neighborRotatedPos[0],3.0)*neighborRotatedPos[1];
            mat(1,1) += std::pow(neighborRotatedPos[0],2.0)*std::pow(neighborRotatedPos[1],2.0);
            mat(1,2) += 0.5*neighborRotatedPos[0]*std::pow(neighborRotatedPos[1],3.0);
            mat(2,0) += 0.25*std::pow(neighborRotatedPos[0],2.0)*std::pow(neighborRotatedPos[1],2.0);
            mat(2,1) += 0.5*neighborRotatedPos[0]*std::pow(neighborRotatedPos[1],3.0);
            mat(2,2) += 0.25*std::pow(neighborRotatedPos[1],4.0);

            b[0] += 0.5*std::pow(neighborRotatedPos[0],2.0)*neighborRotatedPos[2];
            b[1] += neighborRotatedPos[0]*neighborRotatedPos[1]*neighborRotatedPos[2];
            b[2] += 0.5*std::pow(neighborRotatedPos[1],2.0)*neighborRotatedPos[2];
        }
        Eigen::Vector3d ans = mat.bdcSvd(Eigen::ComputeFullU|Eigen::ComputeFullV).solve(b);
        Eigen::EigenSolver<Eigen::Matrix2d> eigensolver;

        Eigen::Matrix2d SecondFundamentalMat;
        SecondFundamentalMat(0,0) = ans[0];
        SecondFundamentalMat(0,1) = ans[1];
        SecondFundamentalMat(1,0) = ans[1];
        SecondFundamentalMat(1,1) = ans[2];

        eigensolver.compute(SecondFundamentalMat);

        #ifdef MY_DEBUG
        std::cout << "2FF matrix of triangle " << SecondFundamentalMat << std::endl;
        #endif 

        Eigen::Vector2d eigenvalues = eigensolver.eigenvalues().real();
        Eigen::Matrix2d eigenvectors= eigensolver.eigenvectors().real();

        Real2 eigenvals;
        eigenvals[0] = eigenvalues[0];
        eigenvals[1] = eigenvalues[1];

        CurvaturePerVertex_[i] = eigenvals;

        // rotate the eigenvectors 
        Real3 eigvec1;
        eigvec1.fill(0);
        Real3 eigvec2;
        eigvec2.fill(0);

        eigvec1[0] = eigenvectors(0,0);
        eigvec1[1] = eigenvectors(1,0);
        eigvec2[0] = eigenvectors(0,1);
        eigvec2[1] = eigenvectors(1,1);

        std::array<Real3, 2> eigenvectorPair;

        Real3 eigvec1Rot = LinAlg3x3::MatrixDotVector(revrotMat, eigvec1);
        Real3 eigvec2Rot = LinAlg3x3::MatrixDotVector(revrotMat, eigvec2);

        // eigenvectors in real space
        principalDir1_[i] = eigvec1Rot;
        principalDir2_[i] = eigvec2Rot;
    }

    for (int i=0;i<CurvaturePerVertex_.size();i++)
    {
        Real avg = 0.0;
        Real gauss = 1.0;
        for (int j=0;j<2;j++)
        {
            avg += CurvaturePerVertex_[i][j]/2.0;
            gauss *= CurvaturePerVertex_[i][j];
        }

        avgCurvaturePerVertex_[i] = avg;
        GaussCurvaturePerVertex_[i] = gauss;
    }
}