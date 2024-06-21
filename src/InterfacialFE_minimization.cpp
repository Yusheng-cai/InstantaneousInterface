#include "InterfacialFE_minimization.h"

namespace MeshRefineStrategyFactory
{
    registry_<InterfacialFE_minimization> registerInterfacialFE_minimization("InterfacialFE_minimization");
}

// InterfacialFE_Solver::InterfacialFE_Solver(std::vector<INT3>& faces) : face_(faces){

// }

// double InterfacialFE_Solver::operator()(const Eigen::VectorXd& xd, Eigen::VectorXd& grad){
//     // convert the information to a mesh
//     std::vector<Real3> verts;
//     for (int i=0;i<3;i++){
//         Real3 p;
//         p[0] = xd[i*3 + 0];
//         p[1] = xd[i*3 + 1];
//         p[2] = xd[i*3 + 2];

//         verts.push_back(p);
//     }

//     Mesh m(verts, face_);

//     MeshTools::CalculateAreaDerivatives(Mesh &m, int &dAdpi)
// }

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

    // boundary terms
    pack_.ReadNumber("boundarymaxstep", ParameterPack::KeyType::Optional, maxBoundaryStep_);
    pack_.ReadNumber("boundarystepsize", ParameterPack::KeyType::Optional,boundarystepsize_);
    pack_.ReadNumber("boundarytolerance", ParameterPack::KeyType::Optional, boundarytol_);
    pack_.ReadNumber("L2maxstep", ParameterPack::KeyType::Optional, maxL2step_);
    pack_.ReadNumber("boundary_optimize_every", ParameterPack::KeyType::Optional, boundary_optimize_every);
    pack_.ReadNumber("dgamma_gamma", ParameterPack::KeyType::Required, dgamma_gamma_);
    pack_.ReadNumber("L2_guess", ParameterPack::KeyType::Optional, L2_);
    pack_.ReadNumber("zstar", ParameterPack::KeyType::Optional, zstar_);
    pack_.ReadNumber("L2_step_size", ParameterPack::KeyType::Optional, L2_stepsize_);
    pack_.ReadNumber("L2tolerance", ParameterPack::KeyType::Optional, L2tol_);
    pack_.ReadNumber("zstar_deviation", ParameterPack::KeyType::Optional, zstar_deviation_);
    pack_.Readbool("useNumerical", ParameterPack::KeyType::Optional, useNumerical_);

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

        // define some variables
        Real max=std::numeric_limits<Real>::lowest();
        Real avg_step = 0;
        int total_verts=0;

        // start updating
        #pragma omp parallel
        {
            Real local_max = std::numeric_limits<Real>::lowest();
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

void InterfacialFE_minimization::refineBoundary(Mesh& m, AFP_shape* shape){
    area_list_.clear(); volume_list_.clear();
    Vnbs_list_.clear(); Anbs_list_.clear();
    std::vector<Real> contact_angle_list, vecArea;
    std::vector<Real3> Normal;

    // first do we refine --> refine updates the mesh
    this->refine(m);

    Real3 Volume_shift={0,0,0};

    if (m.isPeriodic()){
        Volume_shift = -0.5 * m.getBoxLength();
    }

    // outer loop for Lagrange refinement
    int L2_step = 0;
    while (true) {
        Real mean_z;
        Real mean_ca;

        // solve for pi*
        int iteration = 0;
        int cont_ind  = 0;

        Mesh curr_m = m;

        // inner loop for pi, pib, L2 refinement
        while (true){
            // check if we are optimizing mesh
            if ((cont_ind+1) % boundary_optimize_every == 0){
                // then optimize mesh
                MeshTools::CGAL_optimize_Mesh(curr_m,10,60);
                MeshTools::ChangeWindingOrder(curr_m);
            }


                                        // boundary update //
            // calculate area derivative --> dAdr at the boundary and volume derivative dVdr 
            std::vector<Real3> dAdr ,dVdr;
            MeshTools::CalculateAreaDerivatives(curr_m, dAdr);
            MeshTools::CalculateVolumeDerivatives(curr_m, dVdr, Volume_shift);

            // calculate drdu, drdv, boundaryindices, dAnbsdv, dAnbsdu, dVnbsdv, dVnbsdu
            std::vector<int> BoundaryIndices;
            std::vector<Real2> dAnbsduv, dVnbsduv;
            std::vector<Real3> drdu, drdv;
            std::vector<Real> ulist, vlist;
            MeshTools::CalculatedAVnbsdUV(curr_m, shape, BoundaryIndices, ulist, \
                                    vlist, drdu, drdv, dAnbsduv, dVnbsduv, useNumerical_, Volume_shift);

            // access to the vertices in m
            contact_angle_list.clear();

            // set max boundary step to be lower initially = -3.40282e38
            Real max_boundary_step  = std::numeric_limits<Real>::lowest();
            auto& verts             = curr_m.accessvertices();
            int N                   = BoundaryIndices.size();
            Real kk                 = rho_ * (L_ + mu_) / (2*gamma_);
            mean_z                  = 0.0;
            std::cout << "L1 = " << L_ << std::endl;
            std::cout << "L2 = " << L2_ << std::endl;
            std::cout << "k = "  << kk << std::endl;

            // pib refinement
            for (int j=0;j<BoundaryIndices.size();j++){
                // get the actual index of boundary
                int ind       = BoundaryIndices[j];

                Real dAdu     = LinAlg3x3::DotProduct(drdu[j], dAdr[ind]);
                Real dAdv     = LinAlg3x3::DotProduct(drdv[j], dAdr[ind]);

                Real dVdu     = LinAlg3x3::DotProduct(drdu[j], dVdr[ind]);
                Real dVdv     = LinAlg3x3::DotProduct(drdv[j], dVdr[ind]);

                // calculate dAnbsdu and dAnbsdv --> keep drdv the same 
                Real dAnbsdu  = dAnbsduv[j][0];
                Real dAnbsdv  = dAnbsduv[j][1];
                Real dVnbsdu  = dVnbsduv[j][0];
                Real dVnbsdv  = dVnbsduv[j][1];

                // calculate dEdv
                Real dEdv     = dAdv - rho_ * (L_ + mu_) / gamma_* (dVdv + dVnbsdv) + dgamma_gamma_ * dAnbsdv + L2_ * drdv[j][2] / (Real)N;

                // calculate the inverse jacobian
                auto invjac   = shape->InvNumericalJacobian(ulist[j], vlist[j]);
                Eigen::MatrixXd dEduv(2,1);
                dEduv << 0, dEdv;
                auto dEdr     = invjac.transpose() * dEduv;
                Real3 step; step[0]=dEdr(0,0); step[1]=dEdr(1,0); step[2]=dEdr(2,0);

                // we can calculate the contact angle by finding the dgamma_gamma where dEdv is 0
                Real ca = 1.0 / dAnbsdv * (-dAdv + rho_ * (L_ + mu_) / gamma_ * (dVdv + dVnbsdv));
                contact_angle_list.push_back(ca);

                // update vlist 
                verts[ind].position_ = verts[ind].position_ - boundarystepsize_ * step;

                // update the mean z
                mean_z += verts[ind].position_[2];

                if (LinAlg3x3::norm(step) > max_boundary_step){
                    max_boundary_step = LinAlg3x3::norm(step);
                }
            }

            // then do pi refinement
            this->refine(curr_m);

            // calculate area and volume after boundary steps
            vecArea.clear(); Normal.clear();
            a_ = MeshTools::CalculateArea(curr_m, vecArea, Normal);
            V_ = MeshTools::CalculateVolumeDivergenceTheorem(curr_m, vecArea, Normal);
            MeshTools::CalculateAVnbs(curr_m, shape, BoundaryIndices,
                                      ulist, vlist, Anbs_, Vnbs_,10000, useNumerical_, Volume_shift);

            area_list_.push_back(a_);
            volume_list_.push_back(V_);
            Vnbs_list_.push_back(Vnbs_);
            Anbs_list_.push_back(Anbs_);

            // calculate mean z
            mean_z    = mean_z / (Real)N;
            Real var  = Algorithm::calculateVariance(contact_angle_list);
            mean_ca   = Algorithm::calculateMean(contact_angle_list);
            std::cout << "Mean z = " << mean_z << std::endl;
            std::cout << "maxmimum boundary step = " << max_boundary_step << std::endl;
            std::cout << "std of contact angle is " << std::sqrt(var) << std::endl;
            std::cout << "mean of contact angle is " << mean_ca << std::endl; 

            if (iteration > maxBoundaryStep_){
                break;
            }

            // update iterations
            cont_ind++;
            iteration++;

            // check if we deviated too much from zstar
            if (std::abs(mean_z - zstar_) > zstar_deviation_){
                break;
            }

            if (max_boundary_step < boundarytol_){
                break;
            }

        }

        Real L2_step = (mean_z - zstar_);
        L2_list_.push_back(L2_);

        // break the while loop
        if (std::abs(L2_step) < L2tol_ || L2_step > maxL2step_){
            m = curr_m;
            break;
        }
        L2_ = L2_ + L2_stepsize_ * L2_step;

        L2_step++;
    }
}