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
    pack_.ReadNumber("L", ParameterPack::KeyType::Optional, L_);
    pack_.ReadNumber("temperature", ParameterPack::KeyType::Optional, temperature_);
    pack_.ReadNumber("print_every", ParameterPack::KeyType::Optional, print_every);
    pack_.ReadNumber("tolerance", ParameterPack::KeyType::Optional, tol_);
    pack_.ReadNumber("optimize_every", ParameterPack::KeyType::Optional, optimize_every);
    pack_.Readbool("MaxStepCriteria", ParameterPack::KeyType::Optional, MaxStepCriteria);

    // calculate mu and gamma based on temperature
    mu_ = CalculateMu(temperature_);
    gamma_ = CalculateGamma(temperature_);
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

Real InterfacialFE_minimization::CalculateGamma(Real temperature){
    return (30.04 - 0.27477 * (270 - temperature)) * 1e-6 * 1e-18; //kJ/nm2
}

Real InterfacialFE_minimization::CalculateMu(Real temperature){
    return 0.021 * temperature - 5.695;
}


void InterfacialFE_minimization::refine(Mesh& mesh){
    // define mesh
    mesh_ = &mesh;
    mesh_->CalcVertexNormals();

    FE_.clear();

    // first generate necessary things for the mesh
    update_Mesh();

    // calculate cotangent weights 
    for (int i=0;i<maxStep_;i++){
        // calculate the derivatives dAdpi and dVdpi
        MeshTools::CalculateCotangentWeights(*mesh_, neighborIndices_, MapEdgeToFace_, MapEdgeToOpposingVerts_, dAdpi_);
        MeshTools::CalculateVolumeDerivatives(*mesh_, MapVertexToFace_, dVdpi_);
    
        // obtain the vertices
        auto& verts = mesh_->accessvertices();

        Real max=-1e10;
        Real avg_step = 0;
        int total_verts=0;

        // start updating
        #pragma omp parallel
        {
            Real local_max = -1e10;
            Real sum = 0;
            int n_verts= 0 ;
            #pragma omp for
            for (int j=0;j<numVerts_;j++){
                if (! MeshTools::IsBoundary(j, boundaryIndicator_)){
                    // calculate gradient 
                    Real3 gradient = dAdpi_[j] - dVdpi_[j] * rho_*(mu_ + L_) / gamma_;

                    // update position
                    verts[j].position_ = verts[j].position_ - stepsize_ * gradient;

                    // calculate size of step
                    Real step = std::sqrt(LinAlg3x3::DotProduct(gradient, gradient));
                    sum += step;
                    n_verts+=1;

                    if (step > local_max){
                        local_max = step;
                    }
                }
            }
           
            #pragma omp critical
            if (local_max > max){
                max = local_max;
            }

            #pragma omp critical
            avg_step += sum;

            #pragma omp critical
            total_verts += n_verts;
        }

        avg_step /= total_verts;


        // calculate the vertex normals 
        mesh_->CalcVertexNormals();

        if (MaxStepCriteria){
            if (max < tol_){break;}
        }
        else{
            if (avg_step < tol_){break;}
        }

        // print if necessary
        if ((i+1) % print_every == 0){
            std::vector<Real3> Normal;
            std::vector<Real> vecArea;
            Real a = MeshTools::CalculateArea(*mesh_, vecArea, Normal);
            Real V = MeshTools::CalculateVolumeDivergenceTheorem(*mesh_, vecArea, Normal);
            Real E = a - rho_ * (mu_ + L_) / gamma_ * V;
            FE_.push_back(a - rho_ * (mu_ + L_) / gamma_ * V);
            std::cout << "At iteration " << i+1 << " Area = " << a << " " << " Volume = " << V << " Energy = " << E << std::endl;
            std::cout << "max step = " << max << std::endl;
            std::cout << "Avg step = " << avg_step << std::endl;
        }

        if ((i+1) % optimize_every == 0){
            MeshTools::CGAL_optimize_Mesh(*mesh_, 10, 60);

            update_Mesh();
        }
    }

    mesh_->CalcVertexNormals();
}