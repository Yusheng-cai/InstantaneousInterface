#include "CurvatureCurveFit.h"

namespace CurvatureRegistry
{
    registry<CurvatureCurveFit> registerCurvefit("curvefit");
}

CurvatureCurveFit::CurvatureCurveFit(CurvatureInput& input)
:Curvature(input)
{
    input.pack.ReadNumber("neighbors", ParameterPack::KeyType::Optional, NumNeighbors_);

    outputs_.registerOutputFunc("curvature", [this](std::string name) -> void { this -> printCurvature(name);});
}

void CurvatureCurveFit::calculate()
{
    mesh_.findVertexNeighbors();

    const auto& VertexNeighbors_ = mesh_.getNeighborIndices();
    const auto& vertices = mesh_.getvertices();
    std::vector<std::vector<int>> NeighborIndicesNVertex_;
    Graph::getNearbyIndicesNVertexAway(VertexNeighbors_, NumNeighbors_,NeighborIndicesNVertex_);

    CurvaturePerVertex_.clear();
    CurvaturePerVertex_.resize(vertices.size());

    // the reference direction of the normal vector is the z vector
    Real3 referenceDir = {{0,0,1}};

    for (int i=0;i<vertices.size();i++)
    {
        auto& v = vertices[i];
        auto& neighbors = NeighborIndicesNVertex_[i];

        ASSERT((neighbors.size()>=3), "The number of points to be fit must be larger or equal to 3.");

        // How to rotate normal vector to be z vector
        Matrix rotationMat = LinAlg3x3::GetRotationMatrix(v.normals_, referenceDir); 
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
        Real2 eigenvals;
        eigenvals[0] = eigenvalues[0];
        eigenvals[1] = eigenvalues[1];

        CurvaturePerVertex_[i] = eigenvals;
    }
}

void CurvatureCurveFit::printCurvature(std::string name)
{
    std::ofstream ofs_;
    ofs_.open(name);

    ofs_ << "# k1 k2" << "\n";
    for(int i=0;i<CurvaturePerVertex_.size();i++)
    {
        ofs_ << CurvaturePerVertex_[i][0] << " " << CurvaturePerVertex_[i][1] << "\n";
    }
    ofs_.close();
}