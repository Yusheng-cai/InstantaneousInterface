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
    input.pack.ReadString("solver", ParameterPack::KeyType::Optional, solverName_);
    input.pack.Readbool("scale", ParameterPack::KeyType::Optional, scale_);
    input.pack.ReadNumber("k0", ParameterPack::KeyType::Optional, k0_);

    Eigen::initParallel();
    Eigen::setNbThreads(0);
    std::cout << "Eigen is using " << Eigen::nbThreads() << " threads." << std::endl;
}

void MeshCurvatureflow::refine()
{
    mesh_.MapEdgeToFaces();

    // find the boundary vertices
    mesh_.findBoundaryVertices();

    // find map from vertex to the face indices 
    mesh_.MapVertexToFaces();

    // calculate neighbors
    mesh_.findVertexNeighbors();

    if (scale_)
    {
        initialVolume_ = mesh_.calculateVolume();
        std::cout << "Initial volume of the mesh is " << initialVolume_ <<std::endl;
    }

    // find number of vertices 
    const auto& v = mesh_.getvertices();
    newVertices_.insert(newVertices_.end(), v.begin(), v.end());
    TotalArea_.resize(v.size());

    for (int i=0;i<numIterations_;i++)
    {
        std::cout << "Iteration " << i << std::endl;
        if (solverName_ == "explicit")
        {
            refineStep();
        }
        else
        {
            refineImplicitStep();
        }
    }

    mesh_.updateNormals();
}

std::vector<MeshCurvatureflow::Real> MeshCurvatureflow::calculateWeights(int i, std::vector<int>& neighborId)
{
    // obtain the edges that corresponds to this particular vertex 
    auto& edges = mesh_.getEdgeForVertex(i);
    const auto& triangles = mesh_.gettriangles();
    const auto& vertices = mesh_.getvertices();
    int edgesize = edges.size();
    const auto& neighbors = mesh_.getNeighborIndices();

    std::vector<Real> factors;
    Real factor_sum = 0.0;

    // iterate over the edges 
    for (int j=0;j<edgesize;j++)
    {
        // obtain the edge 
        auto& e = edges[j];
        int index1 = e.vertex1_.index;
        int index2 = e.vertex2_.index;

        #ifdef MY_DEBUG
        std::cout << "edge " << j << " index1 = " << index1 << std::endl;
        std::cout << "edge " << j << " index2 = " << index2 << std::endl;
        std::cout << "neighbor edge = " << neighbors[i][j] << std::endl;
        #endif

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
        std::vector<int>& faceIndices = mesh_.getFaceIndicesForEdge(e);

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

            mesh_.getVertexDistance(vertices[index1], vertices[pidx], vec1, vec1sq);
            mesh_.getVertexDistance(vertices[index2], vertices[pidx], vec2, vec2sq);
            mesh_.getVertexDistance(vertices[index1], vertices[index2], vec3, vec3sq);

            Real cos1 = LinAlg3x3::findSinangle(vec1,vec3);
            Real cos2 = LinAlg3x3::findSinangle(vec2,vec3);
            #ifdef MY_DEBUG
            std::cout << "cos1 = " << std::sqrt(1 - cos1*cos1) << std::endl;
            std::cout << "cos2 = " << std::sqrt(1 - cos2*cos2) << std::endl;
            std::cout << "vec1 dist = " << LinAlg3x3::norm(vec1) << std::endl;
            std::cout << "Vec2 dist = " << LinAlg3x3::norm(vec2) << std::endl;
            std::cout << "Vec3 dist = " << LinAlg3x3::norm(vec3) << std::endl;
            #endif

            Real costheta = LinAlg3x3::findCosangle(vec1, vec2);
            Real sintheta = LinAlg3x3::findSinangle(vec1, vec2);

            #ifdef DEBUG
            std::cout << "cosine = " << costheta << std::endl;
            std::cout << "sinetheta = " << sintheta << std::endl;
            std::cout << "Sum of cos2theta and sin2theta is " << costheta*costheta + sintheta*sintheta << std::endl;
            #endif 

            factor += costheta/sintheta;
        }

        factors.push_back(factor);
    }

    return factors;
}

void MeshCurvatureflow::refineStep()
{
    mesh_.CalcTriangleAreaAndFacetNormals();

    // calculate vertex normals but unweighted
    mesh_.CalcVertexNormals();

    const auto& vertexIndicesToface = mesh_.getMapVertexToFace();
    const auto& triangleArea = mesh_.getTriangleArea();
    const auto& vertices = mesh_.getvertices();
    const auto& triangles= mesh_.gettriangles();

    std::fill(TotalArea_.begin(), TotalArea_.end(), 0.0);

    #pragma omp parallel for
    for (int i=0;i<vertexIndicesToface.size();i++)
    {
        for (int j=0;j<vertexIndicesToface[i].size();j++)
        {
            TotalArea_[i] += triangleArea[vertexIndicesToface[i][j]];
        }
    }

    #pragma omp parallel for
    for (int i=0;i<vertices.size();i++)
    {
        // check if the point is on the boundary 
        if (! mesh_.isBoundary(i))
        {
            // obtain the edges that corresponds to this particular vertex 
            auto& edges = mesh_.getEdgeForVertex(i);
            int edgesize = edges.size();

            // initialize the step 
            Real3 step;
            step.fill(0);

            Real factor_sum = 0.0;

            // iterate over the edges 
            for (int j=0;j<edgesize;j++)
            {
                // obtain the edge 
                auto& e = edges[j];
                int index1 = e.vertex1_.index;
                int index2 = e.vertex2_.index;

                #ifdef MY_DEBUG
                std::cout << "edge " << j << " index1 = " << index1 << std::endl;
                std::cout << "edge " << j << " index2 = " << index2 << std::endl;
                #endif

                std::array<int,2> vIndices = {{index1, index2}};
                std::vector<int> otherPoints;

                // obtain the face indices for this particular edge 
                std::vector<int>& faceIndices = mesh_.getFaceIndicesForEdge(e);

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

                Real factor = 0.0;

                for (int k=0;k<otherPoints.size();k++)
                {
                    int pidx = otherPoints[k];
                    Real3 vec1, vec2;
                    Real vec1sq, vec2sq;

                    mesh_.getVertexDistance(vertices[index1], vertices[pidx], vec1, vec1sq);
                    mesh_.getVertexDistance(vertices[pidx], vertices[index2], vec2, vec2sq);
                    

                    Real costheta = LinAlg3x3::findCosangle(vec1, vec2);
                    Real sintheta = LinAlg3x3::findSinangle(vec1, vec2);

                    factor += costheta/sintheta;
                }

                // find the difference of the edge vector
                Real3 diff;
                diff.fill(0);
                for (int k=0;k<3;k++)
                {
                    diff[k] = vertices[index1].position_[k] - vertices[index2].position_[k]; 
                }

                if (e.vertex1_.index == i)
                {
                    for (int k=0;k<3;k++)
                    {
                        diff[k] = -diff[k];
                    }
                }

                for (int k=0;k<3;k++)
                {
                    step[k] += factor * diff[k];
                }

                factor_sum += factor;
            }
            
            #ifdef MY_DEBUG
            std::cout << "Step " << i << " : ";

            for (int j=0;j<3;j++)
            {
                std::cout << step[j] << " ";
            }
            std::cout << "\n";
            #endif

            // multiply step by the factors 
            for (int j=0;j<3;j++)
            {
                newVertices_[i].position_[j] = vertices[i].position_[j] + lambdadt_ * k0_ * vertices[i].normals_[j] + 1.0/(factor_sum) * lambdadt_ * step[j];
            }
        }
    }

    auto& vert = mesh_.accessvertices();
    vert.clear();
    vert.insert(vert.end(), newVertices_.begin(), newVertices_.end());
    mesh_.update();
    std::cout << "Mesh updated." << std::endl;
}

void MeshCurvatureflow::refineImplicitStep()
{
    // obtain the areas of the triangles 
    mesh_.CalcTriangleAreaAndFacetNormals();

    // obtain the normals of the surfaces 
    mesh_.CalcVertexNormalsAreaWeighted();

    auto& vertexIndicesToface = mesh_.getMapVertexToFace();
    const auto& triangleArea = mesh_.getTriangleArea();
    std::fill(TotalArea_.begin(), TotalArea_.end(), 0.0);

    #pragma omp parallel for
    for (int i=0;i<vertexIndicesToface.size();i++)
    {
        for (int j=0;j<vertexIndicesToface[i].size();j++)
        {
            TotalArea_[i] += triangleArea[vertexIndicesToface[i][j]];
        }
    }

    getImplicitMatrix();
    const auto& vertices = mesh_.getvertices();

    // obtain the rhs of the equation to be solved 
    std::vector<Eigen::VectorXf> rhs;
    std::vector<Eigen::VectorXf> xyz;

    rhs.resize(3, Eigen::VectorXf::Zero(vertices.size()));
    xyz.resize(3, Eigen::VectorXf::Zero(vertices.size()));

    #pragma omp parallel for
    for(int i = 0; i < vertices.size(); ++i)
    {
        rhs[0][i] = vertices[i].position_[0];
        rhs[1][i] = vertices[i].position_[1];
        rhs[2][i] = vertices[i].position_[2];
    }

    // solver 
    for (int i=0;i<3;i++)
    {
        xyz[i] = solver_.solveWithGuess(rhs[i],rhs[i]);
    }

    #pragma omp parallel for
    for (int i=0;i<vertices.size();i++)
    {
        for (int j=0;j<3;j++)
        {
            newVertices_[i].position_[j] = xyz[j][i];
        }
    }

    auto& vert = mesh_.accessvertices();
    vert.clear();
    vert.insert(vert.end(), newVertices_.begin(), newVertices_.end());

    // calculate the volume if needed 
    if (scale_)
    {
        Real vol = mesh_.calculateVolume();
        std::cout << "Volume = " << vol << std::endl;
        Real scale = std::pow(initialVolume_/vol, 1.0/3.0);

        // this already performs update 
        mesh_.scaleVertices(scale);
    }
    else
    {
        mesh_.update();
    }
}

void MeshCurvatureflow::getImplicitMatrix()
{
    const auto& vertices = mesh_.getvertices();
    L_.setZero();
    L_.resize(vertices.size(), vertices.size());

    const auto& neighborIndices = mesh_.getNeighborIndices();
    Real epsilon=1e-5;

    triplets_.clear();
    triplets_buffer_.set_master_object(triplets_);

    #pragma omp parallel
    {
        auto& trip_omp = triplets_buffer_.access_buffer_by_id();
        trip_omp.clear();

        #pragma omp for
        for (int i=0;i<vertices.size();i++)
        {
            if (! mesh_.isBoundary(i))
            {
                int numneighbors = neighborIndices[i].size();
                std::vector<int> neighborId;

                // get the weights of the neighbor indices 
                std::vector<Real> weights = calculateWeights(i, neighborId);
                ASSERT((weights.size() == numneighbors), "The number of weights provided for a vertex does not agree with the number of neighbors.");
                Real w = 0.0;

                // Skip the calculation is Total area is 0 --> less than some threshold --> 1e-5
                if (TotalArea_[i] > epsilon)
                {
                    Real factor = 1.0/(4.0 * TotalArea_[i]);

                    for (int j=0;j<numneighbors;j++)
                    {
                        w += weights[j];
                    }

                    for (int j=0;j<numneighbors;j++)
                    {
                        trip_omp.push_back(triplet(i, neighborId[j], factor * weights[j]));
                    }

                    trip_omp.push_back(triplet(i,i, -factor * w));
                }
            }
        }
    }

    int size = triplets_.size();
    for (auto it = triplets_buffer_.beginworker(); it < triplets_buffer_.endworker(); it++)
    {
        size += it->size();
    }
    triplets_.reserve(size);
    for (auto it = triplets_buffer_.beginworker(); it < triplets_buffer_.endworker(); it++)
    {
        triplets_.insert(triplets_.end(), it -> begin(), it -> end());
    }

    L_.setFromTriplets(triplets_.begin(), triplets_.end());

    Eigen::SparseMatrix<Real> I = Eigen::MatrixXf::Identity(vertices.size(), vertices.size()).sparseView();

    L_ = I - lambdadt_*L_;

    // Solve for x, y, z
    // solver_.analyzePattern(L_);
    // solver_.factorize(L_);
    solver_.compute(L_);
}