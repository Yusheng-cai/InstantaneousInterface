#include "CurvatureJetFit.h"

namespace CurvatureRegistry
{
    registry<CurvatureJetFit> registerJet("jetfit");
}

CurvatureJetFit::CurvatureJetFit(CurvatureInput& input)
:Curvature(input)
{
    input.pack.ReadNumber("neighbors", ParameterPack::KeyType::Required, numneighbors_);
    input.pack.ReadNumber("degree", ParameterPack::KeyType::Optional, degree_);
    input.pack.ReadNumber("MongeCoefficient", ParameterPack::KeyType::Optional, MongeCoefficient_);

    outputs_.registerOutputFunc("neighborIndices", [this](std::string name) -> void {this -> printNeighbors(name);});
    outputs_.registerOutputFunc("coefficients", [this](std::string name) -> void {this -> printCoefficientPerVertex(name);});
    outputs_.registerOutputFunc("PCAeigenvector", [this](std::string name) -> void {this -> printPCAeigenvector(name);});

    numpoints_ = (degree_ + 1)*(degree_+2)/2;
}

void CurvatureJetFit::printPCAeigenvector(std::string name)
{
    std::ofstream ofs_;
    ofs_.open(name);

    for (int i=0;i<PCAeigenvector_.size();i++)
    {
        for (int j=0;j<3;j++)
        {
            for (int k=0;k<3;k++)
            {
                ofs_ << PCAeigenvector_[i][j][k] << " ";
            }
        }
        ofs_ << "\n";
    }
    ofs_.close();
}

void CurvatureJetFit::calculate()
{
    std::cout << "Calculating Jet fitting." << std::endl;
    mesh_.findVertexNeighbors();

    const auto& VertexNeighbors_ = mesh_.getNeighborIndices();
    const auto& vertices = mesh_.getvertices();
    Graph::getNearbyIndicesNVertexAway(VertexNeighbors_,numneighbors_,NeighborIndicesNVertex_);
    CurvaturePerVertex_.resize(vertices.size());
    coefficientsPerVertex_.resize(vertices.size());
    PCAeigenvector_.resize(vertices.size());

    for (int i=0;i<NeighborIndicesNVertex_.size();i++)
    {
        std::cout << i << " has " << NeighborIndicesNVertex_[i].size() << " neighbors." << std::endl;
        std::vector<Dpoint> vec_;
        auto& thisVertex = vertices[i];
        auto& thisPos = vertices[i].position_;
        auto& thisNormal = thisVertex.normals_;
        Dpoint point(thisPos[0], thisPos[1], thisPos[2]);
        vec_.push_back(point);

        for (int j=0;j<NeighborIndicesNVertex_[i].size();j++)
        {
            int neighborIndex = NeighborIndicesNVertex_[i][j];
            auto& vertpos = vertices[neighborIndex].position_;

            Dpoint point(vertpos[0], vertpos[1], vertpos[2]);
            vec_.push_back(point);
        }

        mform_ = jetfitter_(vec_.begin(), vec_.end(), degree_, MongeCoefficient_);
        for (int j=0;j<3;j++)
        {
            for (int k=0;k<3;k++)
            {
                PCAeigenvector_[i][j][k] = jetfitter_.pca_basis(j).second[k];
            }
        }

        DVector norm(thisNormal[0], thisNormal[1], thisNormal[2]);

        mform_.comply_wrt_given_normal(norm);
        #ifdef MY_DEBUG
        std::cout << "Origin = " << mform_.origin() << std::endl;
        std::cout << "Position1 = " << vec_[0] << std::endl;
        std::cout << "vertex position = " << vertices[i].position_[0] << " " << vertices[i].position_[1] << " " << vertices[i].position_[2] << std::endl;
        #endif 

        std::vector<Real> coeff_;
        for (int i=0;i<mform_.coefficients().size();i++)
        {
            coeff_.push_back(mform_.coefficients()[i]);
        }


        coefficientsPerVertex_[i] = coeff_;
        CurvaturePerVertex_[i][0] = mform_.coefficients()[0];
        CurvaturePerVertex_[i][1] = mform_.coefficients()[1];
    }
}


void CurvatureJetFit::printCurvature(std::string name)
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

void CurvatureJetFit::printNeighbors(std::string name)
{
    std::ofstream ofs_;
    ofs_.open(name);

    for (int i=0;i<NeighborIndicesNVertex_.size();i++)
    {
        int neighborsize = NeighborIndicesNVertex_[i].size();
        // print self first 
        ofs_ << i << " ";

        for (int j=0;j<neighborsize;j++)
        {
            ofs_ << NeighborIndicesNVertex_[i][j] << " ";
        }
        ofs_ << "\n";
    }
    ofs_.close();
}

void CurvatureJetFit::printCoefficientPerVertex(std::string name)
{
    std::ofstream ofs_;
    ofs_.open(name);

    for (int i=0;i<coefficientsPerVertex_.size();i++)
    {
        int sizePerVertex = coefficientsPerVertex_[i].size();
        for (int j=0;j<sizePerVertex;j++)
        {
            ofs_ << coefficientsPerVertex_[i][j] << " ";
        }
        ofs_ << "\n";
    }
    ofs_.close();
}