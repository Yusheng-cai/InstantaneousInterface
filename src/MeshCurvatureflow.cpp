#include "MeshCurvatureflow.h"

namespace MeshRefineStrategyFactory
{
    registry_<MeshCurvatureflow> registercurveflow("curvatureflow");
}

MeshCurvatureflow::MeshCurvatureflow(MeshRefineStrategyInput& input)
: MeshRefineStrategy(input)
{
    input.pack.ReadNumber("iterations", ParameterPack::KeyType::Required, numIterations_);
    input.pack.ReadNumber("lambdadt", ParameterPack::KeyType::Optional, lambdadt_);
    input.pack.Readbool("scale", ParameterPack::KeyType::Optional, scale_);
    input.pack.Readbool("virtualSite", ParameterPack::KeyType::Optional, virtualSite_);

    Eigen::initParallel();
    Eigen::setNbThreads(0);
    std::cout << "Eigen is using " << Eigen::nbThreads() << " threads." << std::endl;
}

void MeshCurvatureflow::refine(Mesh& mesh)
{
    mesh_ = &mesh;

    MeshTools::MapEdgeToFace(*mesh_, MapEdgeToFace_, MapVertexToEdge_);
    MeshTools::CalculateBoundaryVertices(*mesh_, MapEdgeToFace_, boundaryIndicator_);
    MeshTools::CalculateVertexNeighbors(*mesh_, neighborIndices_);
    MeshTools::MapVerticesToFaces(*mesh_, MapVertexToFace_);

    if (scale_)
    {
        initialVolume_ = mesh_->calculateVolume();
        std::cout << "Initial volume of the mesh is " << initialVolume_ <<std::endl;
    }

    // find number of vertices 
    const auto& v = mesh_->getvertices();
    numVerts_ = v.size();
    newVertices_.resize(numVerts_);
    TotalArea_.resize(numVerts_);

    // start refining 
    for (int i=0;i<numIterations_;i++)
    {
        std::cout << "Iteration " << i << std::endl;
        refineImplicitStep();
    }
}

std::vector<MeshCurvatureflow::Real> MeshCurvatureflow::calculateWeights(int i, std::vector<int>& neighborId, Real3& Lfactor, bool& flag)
{
    // obtain triangles and vertices from mesh
    const auto& triangles = mesh_->gettriangles();
    const auto& vertices  = mesh_->getvertices();

    // obtain the edges that corresponds to this particular vertex 
    auto edgesIndices     = MapVertexToEdge_[i];
    int edgesize = edgesIndices.size();
    ASSERT((edgesize > 0), "Vertex " << i << " corresponds to 0 edges.");

    // factor at which we detect INF
    Real epsilon=1e-8;

    std::vector<Real> factors;
    Real factor_sum = 0.0;

    flag = true;

    // iterate over the edges 
    for (int j=0;j<edgesize;j++)
    {
        // obtain the edge 
        auto& e = edgesIndices[j];
        int index1 = e[0];
        int index2 = e[1];

        if (index1 != i)
        {
            neighborId.push_back(index1);
        }
        else
        {
            neighborId.push_back(index2);
        }

        std::array<int,2> vIndices = {{index1, index2}};

        // the other 2 points out of the 4 points that forms 2 adjacent triangles
        std::vector<int> otherPoints;

        // obtain the face indices for this particular edge 
        auto it = MapEdgeToFace_.find(e);
        ASSERT((it != MapEdgeToFace_.end()), "The edge " << e[0] << " " << e[1] << " is not found.");
        std::vector<int> faceIndices  = it -> second;
        ASSERT((faceIndices.size() == 2), "We are only smoothing the nonboundary points so number of faces per vertex must be 2 while it is " << faceIndices.size());

        for (int k=0;k<faceIndices.size();k++)
        {
            auto& faceidx = triangles[faceIndices[k]].triangleindices_;
            for (int m=0;m<3;m++)
            {
                auto found = (std::find(vIndices.begin(), vIndices.end(), faceidx[m])) != vIndices.end();

                if (! found)
                {
                    otherPoints.push_back(faceidx[m]);
                }
            }
        }

        ASSERT((otherPoints.size() == 2), "There must be 2 other points that forms a pair triangle along with the edge while it is " << otherPoints.size());
        ASSERT((otherPoints[0] != otherPoints[1]), "The 2 triangle indices must be unique.");

        Real factor = 0.0;
        for (int k=0;k<otherPoints.size();k++)
        {
            int pidx = otherPoints[k];
            Real3 vec1, vec2, vec3;
            Real vec1sq, vec2sq, vec3sq;

            mesh_->getVertexDistance(vertices[index1], vertices[pidx], vec1, vec1sq);
            mesh_->getVertexDistance(vertices[index2], vertices[pidx], vec2, vec2sq);
            Real costheta = LinAlg3x3::findCosangle(vec1, vec2);
            Real num = 1 - costheta*costheta;
            if (num < epsilon)
            {
                flag = false;
            }
            Real sintheta = std::sqrt(1 - costheta*costheta);
            
            factor += costheta/sintheta;
        }

        if (mesh_->isPeriodic())
        {
            Real3 ShiftVec;
            mesh_ -> CalculateShift(vertices[neighborId[j]].position_, vertices[i].position_,ShiftVec);

            for (int k=0;k<3;k++)
            {
                Lfactor[k] += factor * ShiftVec[k];
            }
        }

        factors.push_back(factor);
    }

    return factors;
}

void MeshCurvatureflow::refineImplicitStep()
{
    // obtain the normals of the surfaces 
    mesh_->CalcVertexNormals();

    // Get sum of triangles for each of the vertices  
    std::vector<Real3> normals;
    MeshTools::CalculateTriangleAreasAndFaceNormals(*mesh_, TriangleAreas_, normals);
    std::fill(TotalArea_.begin(), TotalArea_.end(), 0.0);

    #pragma omp parallel for
    for (int i=0;i<MapVertexToFace_.size();i++)
    {
        for (int j=0;j<MapVertexToFace_[i].size();j++)
        {
            TotalArea_[i] += TriangleAreas_[MapVertexToFace_[i][j]];
        }
    }

    // obtain the implicit matrices
    getImplicitMatrix();

    // obtain the rhs of the equation to be solved 
    std::vector<Eigen::VectorXf> rhs;
    std::vector<Eigen::VectorXf> xyz;
    rhs.resize(3, Eigen::VectorXf::Zero(numVerts_));
    xyz.resize(3, Eigen::VectorXf::Zero(numVerts_));

    // populate the solved matrices 
    const auto& vertices = mesh_->getvertices();
    #pragma omp parallel for
    for(int i = 0; i < numVerts_; ++i)
    {
        // if mesh is not periodic and we are not doing virtual site
        // then we must scale the rhs by Lfactors_ as in (xj - xi + L)
        rhs[0][i] = lambdadt_ * Lfactors_[i][0] + vertices[i].position_[0];
        rhs[1][i] = lambdadt_ * Lfactors_[i][1] + vertices[i].position_[1];
        rhs[2][i] = lambdadt_ * Lfactors_[i][2] + vertices[i].position_[2];
    }

    // solve the sparse system of equations
    Sparse_LU sparsesolver;
    sparsesolver.compute(L_);
    ASSERT((sparsesolver.info() == Eigen::Success), "Compute step failed.");
    for (int i=0;i<3;i++)
    {
        xyz[i] = sparsesolver.solve(rhs[i]);
        ASSERT((sparsesolver.info() == Eigen::Success), "The solver failed.");
    }

    // copy over to new vertices data 
    #pragma omp parallel for
    for (int i=0;i<numVerts_;i++)
    {
        for (int j=0;j<3;j++)
        {
            newVertices_[i].position_[j] = xyz[j][i];
        }
    }

    // update the mesh vertices 
    auto& vert = mesh_->accessvertices();
    vert.clear();
    vert.insert(vert.end(), newVertices_.begin(), newVertices_.end());

    // calculate the volume if needed 
    if (scale_)
    {
        Real vol = mesh_->calculateVolume();
        Real scale = std::pow(initialVolume_/vol, 1.0/3.0);

        // this already performs update 
        mesh_->scaleVertices(scale);
    }
    else
    {
        mesh_->update();
    }
}

void MeshCurvatureflow::getImplicitMatrix()
{
    // What to do for periodic Meshes?
    /*
        We can make virtual vertices that updates every step.

        For each vertex i, we can check its neighbors, find their distances. If the distance is larger
        than box/2 as usually defined for PBC, we can then make a virtual vertex for vertex i. 
    */
    const auto& vertices = mesh_->getvertices();

    Lfactors_.clear();
    Lfactors_.resize(vertices.size(), {{0,0,0}});
    triplets_.clear();

    #pragma omp parallel
    {
        std::vector<triplet> triplet_local;

        #pragma omp for
        for (int i=0;i<vertices.size();i++)
        {
            if (! MeshTools::IsBoundary(i, boundaryIndicator_))
            {
                int numneighbors = neighborIndices_[i].size();
                std::vector<int> neighborId;

                // get the weights of the neighbor indices 
                Real3 Lfactor = {};
                bool flag = true;
                std::vector<Real> weights = calculateWeights(i, neighborId, Lfactor, flag);
                ASSERT((weights.size() == numneighbors), "The number of weights provided for a vertex does not agree with the number of neighbors.");
                Real w = 0.0;

                // Skip the calculation is Total area is 0 --> less than some threshold --> 1e-5
                if (flag)
                {
                    Real factor = 1.0/(4.0 * TotalArea_[i]);

                    for (int j=0;j<numneighbors;j++)
                    {
                        w += weights[j];
                    }

                    for (int j=0;j<numneighbors;j++)
                    {
                        triplet_local.push_back(triplet(i, neighborId[j], factor * weights[j]));
                    }

                    for (int j=0;j<3;j++)
                    {
                        Lfactors_[i][j] = Lfactor[j] * factor;
                    }

                    triplet_local.push_back(triplet(i,i, -factor * w));
                }
            }
        }

        #pragma omp critical
        {
            triplets_.insert(triplets_.end(), triplet_local.begin(), triplet_local.end());
        }
    }

    // set the L matrix 
    L_.setZero();
    L_.resize(numVerts_, numVerts_);
    L_.setFromTriplets(triplets_.begin(), triplets_.end());
    Eigen::SparseMatrix<Real> I = Eigen::MatrixXf::Identity(numVerts_, numVerts_).sparseView();
    L_ = I - lambdadt_*L_;
}