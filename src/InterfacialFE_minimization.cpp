#include "InterfacialFE_minimization.h"

namespace MeshRefineStrategyFactory
{
    registry_<InterfacialFE_minimization> registerInterfacialFE_minimization("InterfacialFE_minimization");
}

InterfacialFE_minimization::InterfacialFE_minimization(MeshRefineStrategyInput& input)
    : MeshRefineStrategy(input)
{
    pack_.ReadNumber("maxstep", ParameterPack::KeyType::Optional, maxStep_);
    pack_.ReadNumber("stepsize", ParameterPack::KeyType::Optional, stepsize_);
    pack_.ReadNumber("k0", ParameterPack::KeyType::Required, k0_);
}

void InterfacialFE_minimization::update_Mesh(){
    // obtain map from edge index {minIndex, maxIndex} to the face index 
    // obtain map from vertex index {index} to the Edge Index {minIndex, maxIndex}
    MeshTools::MapEdgeToFace(*mesh_, MapEdgeToFace_, MapVertexToEdge_);

    // Calculate the boundary vertices --> vertices which that has an edge shared by only 1 face
    MeshTools::CalculateBoundaryVertices(*mesh_, MapEdgeToFace_, boundaryIndicator_);

    // Calculate the vertex neighbors
    MeshTools::CalculateVertexNeighbors(*mesh_, neighborIndices_);

    // Map from vertex indices to face indices 
    MeshTools::MapVerticesToFaces(*mesh_, MapVertexToFace_);

    // Map from edges to their opposing vertex indices 
    MeshTools::MapEdgeToOpposingVertices(*mesh_, MapEdgeToFace_, MapEdgeToOpposingVerts_);

    // find number of vertices 
    const auto& v = mesh_->getvertices();
    numVerts_ = v.size();
    newVertices_.clear();newVertices_.resize(numVerts_);
    TotalArea_.clear();TotalArea_.resize(numVerts_);
}

void InterfacialFE_minimization::refine(Mesh& mesh){
    // define mesh
    mesh_ = &mesh;
    mesh_->CalcVertexNormals();
    
    // obtain the vertices
    auto& verts = mesh_->accessvertices();

    // first generate necessary things for the mesh
    update_Mesh();

    // calculate cotangent weights 
    for (int i=0;i<maxStep_;i++){
        // calculate the derivatives dAdpi and dVdpi
        MeshTools::CalculateCotangentWeights(*mesh_, neighborIndices_, boundaryIndicator_, MapEdgeToFace_, MapEdgeToOpposingVerts_, dAdpi_);
        MeshTools::CalculateVolumeDerivatives(*mesh_, MapVertexToFace_, dVdpi_);


        // start updating
        #pragma omp parallel for
        for (int i=0;i<numVerts_;i++){
            if (! MeshTools::IsBoundary(i, boundaryIndicator_)){
                verts[i].position_ = verts[i].position_ - stepsize_ * (dAdpi_[i] - dVdpi_[i] * 2.0 * k0_);
            }
        }

        // calcuilate vertex normal again 
        mesh_->CalcVertexNormals();

        // compare vertex normal with dA --> dA = 2 kappa * N
        // for (int i=0;i<dAdpi_.size();i++){
        //     if (! MeshTools::IsBoundary(i,boundaryIndicator_)){
        //         Real3 dA = dAdpi_[i];

        //         Real product = LinAlg3x3::DotProduct(dA, dA);

        //         if (product > 1e-8){
        //             LinAlg3x3::normalize(dA);
        //             std::cout << "Calculated normal = " << dA << " , actual normal = " << mesh_->getvertices()[i].normals_ << std::endl; 
        //         }
        //     }
        // }

        std::vector<Real3> Normal;
        std::vector<Real> vecArea;
        Real a = MeshTools::CalculateArea(*mesh_, vecArea, Normal);
        Real V = MeshTools::CalculateVolumeDivergenceThoerem(*mesh_, vecArea, Normal);
        std::cout << "Area = " << a << std::endl;
        std::cout << "Volume = " << V << std::endl;
    }
}