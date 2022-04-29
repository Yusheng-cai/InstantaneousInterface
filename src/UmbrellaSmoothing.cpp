#include "UmbrellaSmoothing.h"

namespace MeshRefineStrategyFactory
{
    registry_<UmbrellaSmoothing> registerUmbrella("umbrella");
}

UmbrellaSmoothing::UmbrellaSmoothing(MeshRefineStrategyInput& input)
:MeshRefineStrategy(input)
{
    input.pack.ReadNumber("iterations", ParameterPack::KeyType::Required, numIterations_);
    input.pack.ReadNumber("lambdadt", ParameterPack::KeyType::Optional, lambdadt_);
    input.pack.Readbool("scale", ParameterPack::KeyType::Optional, scale_);
}

void UmbrellaSmoothing::refine(Mesh& mesh)
{
    mesh_ = &mesh;

    mesh_->findVertexNeighbors();
    mesh_->findBoundaryVertices();
    const auto& vertices = mesh_->getvertices();

    // fill the old vertices with the original vertices
    oldVertices_.insert(oldVertices_.end(), vertices.begin(), vertices.end());

    // resize the new vertices 
    newVertices_.insert(newVertices_.end(), vertices.begin(), vertices.end());

    // find volume if needed 
    if ( scale_)
    {
        initialVolume_ = mesh_->calculateVolume();
        std::cout << "Initial volume = " << initialVolume_ << std::endl;
    }

    prepareImplicitMatrix();

    for (int i=0;i<numIterations_;i++)
    {
        std::cout << "Step " << i << std::endl;
        refineStepImplicit();
    }

    auto& vert = mesh_->accessvertices();
    vert.clear();
    vert.insert(vert.end(), newVertices_.begin(), newVertices_.end());
    mesh_->updateNormals();
}

void UmbrellaSmoothing::refineStepImplicit()
{
    const auto& vertices = mesh_->getvertices();

    // obtain the rhs of the equation to be solved 
    std::vector<Eigen::VectorXf> rhs;
    std::vector<Eigen::VectorXf> xyz;

    rhs.resize(3, Eigen::VectorXf::Zero(vertices.size()));
    xyz.resize(3, Eigen::VectorXf::Zero(vertices.size()));

    for(int i = 0; i < vertices.size(); ++i)
    {
        rhs[0][i] = vertices[i].position_[0];
        rhs[1][i] = vertices[i].position_[1];
        rhs[2][i] = vertices[i].position_[2];
    }

    for (int i=0;i<3;i++)
    {
        xyz[i] = solver_.solveWithGuess(rhs[i],rhs[i]);
    }

    for (int i=0;i<vertices.size();i++)
    {
        for (int j=0;j<3;j++)
        {
            newVertices_[i].position_[j] = xyz[j][i];
        }
    }

    auto& vert = mesh_->accessvertices();
    vert.clear();
    vert.insert(vert.end(), newVertices_.begin(), newVertices_.end());

    if (scale_)
    {
        Real vol = mesh_->calculateVolume();
        Real scale = std::pow(initialVolume_/vol, 1.0/3.0);
        mesh_->scaleVertices(scale);
    }
    else
    {
        mesh_->update();
    }
}


void UmbrellaSmoothing::prepareImplicitMatrix()
{
    const auto& vertices = mesh_->getvertices();
    L_.resize(vertices.size(), vertices.size());

    const auto& neighborIndices = mesh_->getNeighborIndices();

    triplets_.clear();
    // 10 is a guess here 
    triplets_.reserve(vertices.size()*10);

    for (int i=0;i<vertices.size();i++)
    {
        if (! mesh_->isBoundary(i))
        {
            int numneighbors = neighborIndices[i].size();
            Real factor = 1.0/numneighbors;

            for (int j=0;j<numneighbors;j++)
            {
                triplets_.push_back(triplet(i, neighborIndices[i][j], factor));
            }

            triplets_.push_back(triplet(i,i, -1));
        }
    }
    L_.setFromTriplets(triplets_.begin(), triplets_.end());

    Eigen::SparseMatrix<Real> I = Eigen::MatrixXf::Identity(vertices.size(), vertices.size()).sparseView();

    L_ = I - lambdadt_*L_;

    // Solve for x, y, z
    // solver_.analyzePattern(L_);
    // solver_.factorize(L_);
    solver_.compute(L_);
}