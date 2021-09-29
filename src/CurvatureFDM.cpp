#include "CurvatureFDM.h"

namespace CurvatureRegistry
{
    registry<CurvatureFDM> registerFDM("FDM");
}

CurvatureFDM::CurvatureFDM(CurvatureInput& input)
:Curvature(input)
{
    input.pack.ReadString("mean", ParameterPack::KeyType::Optional, MeanMethod_);
    ASSERT((MeanMethod_ == "arithmetic" || MeanMethod_ == "geometric"), "The method for calculating mean must be either arithmetic or geometric.");

    outputs_.registerOutputFunc("curvature", [this](std::string name) -> void { this -> printCurvature(name);});
}

void CurvatureFDM::calculate()
{
    mesh_.findVertexNeighbors();

    if (MeanMethod_ == "arithmetic")
    {
        calculateArithmeticMean();
    }
    else
    { 
        calculateGeometricMean();
    }
}

void CurvatureFDM::printCurvature(std::string name)
{
    std::ofstream ofs_;
    ofs_.open(name);

    if (MeanMethod_ == "arithmetic")
    {
        ofs_ << "# Mean curvature" << "\n";
    }
    else
    {
        ofs_ << "# Gaussian Curvature" << "\n";
    }
    for (int i=0;i<curvature_.size();i++)
    {
        ofs_ << curvature_[i];
        ofs_ << "\n"; 
    }
    ofs_.close();
}

void CurvatureFDM::calculateArithmeticMean()
{   
    curvature_.clear();
    curvature_.resize(mesh_.getNumVertices());
    std::fill(curvature_.begin(), curvature_.end(),0.0);

    const auto& vertices_ = mesh_.getvertices();
    const auto& neighbor_indices_ = mesh_.getNeighborIndices();
    const auto& numNeighbors = mesh_.getNumNeighbors();

    // Now we calculate the curvature
    for (int i=0;i<vertices_.size();i++)
    {
        auto& v1 = vertices_[i];
        auto& neighbors = neighbor_indices_[i];

        Real neighborAreaTot = 0.0;;

        for (int j=0;j<neighbors.size();j++)
        { 
            auto& v2 = vertices_[neighbors[j]];

            Real3 diff;
            Real3 diffn;
            Real sq_diff=0.0;
            Real dot_product = 0.0;
            Real curve_Val;
            diff.fill(0);
            diffn.fill(0);
            
            for(int m=0;m<3;m++)
            {
                diff[m] = v2.position_[m] - v1.position_[m];
                diffn[m] = v2.normals_[m] - v1.normals_[m];
            }

            Real diffnnorm = LinAlg3x3::norm(diffn); 

            for (int m=0;m<3;m++)
            {
                sq_diff += diff[m]*diff[m];
            }

            for (int m=0;m<3;m++)
            {
                dot_product += diff[m]*diffn[m];
            }

            curve_Val = dot_product/sq_diff;
            curvature_[i]  += std::abs(curve_Val);
        }
    }

    ASSERT((numNeighbors.size() == curvature_.size()), "The size for number of neighors is wrong.");

    for (int i=0;i<curvature_.size();i++)
    { 
        curvature_[i] = curvature_[i]/numNeighbors[i];
    } 

}

void CurvatureFDM::calculateGeometricMean()
{
    curvature_.clear();
    curvature_.resize(mesh_.getNumVertices());
    std::fill(curvature_.begin(), curvature_.end(),1.0);

    const auto& vertices_ = mesh_.getvertices();
    const auto& neighbor_indices_ = mesh_.getNeighborIndices();
    const auto& numNeighbors = mesh_.getNumNeighbors();

    // Now we calculate the curvature
    for (int i=0;i<vertices_.size();i++)
    {
        auto& v1 = vertices_[i];
        auto& neighbors = neighbor_indices_[i];

        Real neighborAreaTot = 0.0;;

        for (int j=0;j<neighbors.size();j++)
        { 
            auto& v2 = vertices_[neighbors[j]];

            Real3 diff;
            Real3 diffn;
            Real sq_diff=0.0;
            Real dot_product = 0.0;
            Real curve_Val;
            diff.fill(0);
            diffn.fill(0);
            
            for(int m=0;m<3;m++)
            {
                diff[m] = v2.position_[m] - v1.position_[m];
                diffn[m] = v2.normals_[m] - v1.normals_[m];
            }

            Real diffnnorm = LinAlg3x3::norm(diffn); 

            for (int m=0;m<3;m++)
            {
                sq_diff += diff[m]*diff[m];
            }

            for (int m=0;m<3;m++)
            {
                dot_product += diff[m]*diffn[m];
            }

            curve_Val = dot_product/sq_diff;
            curvature_[i]  *= std::abs(curve_Val);
        }
    }

    ASSERT((numNeighbors.size() == curvature_.size()), "The size for number of neighors is wrong.");

    for (int i=0;i<curvature_.size();i++)
    { 
        curvature_[i] = std::pow(curvature_[i], 1.0/numNeighbors[i]);
    } 
}