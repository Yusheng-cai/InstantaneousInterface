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

std::vector<MeshCurvatureflow::Real> MeshCurvatureflow::calculateWeights(int i, std::vector<int>& neighborId, Real3& Lfactor)
{
    // obtain the edges that corresponds to this particular vertex 
    auto& edgesIndices = mesh_.getEdgeIndexForVertex(i);
    const auto& triangles = mesh_.gettriangles();
    const auto& vertices = mesh_.getvertices();
    int edgesize = edgesIndices.size();

    std::vector<Real> factors;
    Real factor_sum = 0.0;

    // iterate over the edges 
    for (int j=0;j<edgesize;j++)
    {
        // obtain the edge 
        auto& e = edgesIndices[j];
        int index1 = e[0];
        int index2 = e[1];

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
        std::vector<int>& faceIndices = mesh_.getFaceIndicesForEdgeIndex(e);

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

            #ifdef MY_DEBUG
            mesh_.getVertexDistance(vertices[index1], vertices[index2], vec3, vec3sq);
            Real cos1 = LinAlg3x3::findSinangle(vec1,vec3);
            Real cos2 = LinAlg3x3::findSinangle(vec2,vec3);
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

            #ifdef DEBUG
            if (std::abs(costheta/sintheta) > 10)
            {
                std::cout << "The three triangles are : " << "\n";
                std::cout << index1 << ". " << vertices[index1].position_[0] << " " << vertices[index1].position_[1] << " " << vertices[index1].position_[2] << "\n";
                std::cout << pidx << ". " << vertices[pidx].position_[0] << " " << vertices[pidx].position_[1] << " " << vertices[pidx].position_[2] << "\n";
                std::cout << index2 << ". " << vertices[index2].position_[0] << " " << vertices[index2].position_[1] << " " << vertices[index2].position_[2] << "\n";
                std::cout << "vec1 = " << vec1[0] << " " << vec1[1] << " " << vec1[2] << "\n";
                std::cout << "vec2 = " << vec2[0] << " " << vec2[1] << " " << vec2[2] << "\n";
                std::cout << "Costheta = " << costheta << "\n";
                std::cout << "Sintheta = " << sintheta << "\n";
                Real3 crossproduct = LinAlg3x3::CrossProduct(vec1, vec2);
                Real area = LinAlg3x3::norm(crossproduct);
                std::cout << "Area = " << area << "\n";
            }
            #endif
        }

        if (mesh_.isPeriodic())
        {
            auto boxLength = mesh_.getBoxLength();

            for (int k=0;k<3;k++)
            {
                Real diff = vertices[neighborId[j]].position_[k] - vertices[i].position_[k];
                Real Lj = 0.0;

                if (diff > boxLength[k] * 0.5)
                {
                    Lj = - boxLength[k];
                }

                if (diff < -boxLength[k] * 0.5)
                {
                    Lj = boxLength[k];
                }


                Lfactor[k] += factor * Lj;
            }
        }

        factors.push_back(factor);
    }

    return factors;
}

void MeshCurvatureflow::refineStep()
{
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
            auto& edgeIndices = mesh_.getEdgeIndexForVertex(i);
            int edgesize = edgeIndices.size();

            // initialize the step 
            Real3 step;
            step.fill(0);

            Real factor_sum = 0.0;

            // iterate over the edges 
            for (int j=0;j<edgesize;j++)
            {
                // obtain the edge 
                auto& vIndices = edgeIndices[j];
                int idx1 = vIndices[0];
                int idx2 = vIndices[1];
                std::vector<int> otherPoints;

                #ifdef MY_DEBUG
                std::cout << "edge " << j << " index1 = " << vIndices[0] << std::endl;
                std::cout << "edge " << j << " index2 = " << vIndices[1] << std::endl;
                #endif

                // obtain the face indices for this particular edge 
                std::vector<int>& faceIndices = mesh_.getFaceIndicesForEdgeIndex(vIndices);

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

                    mesh_.getVertexDistance(vertices[idx1], vertices[pidx], vec1, vec1sq);
                    mesh_.getVertexDistance(vertices[idx2], vertices[pidx], vec2, vec2sq);
                    

                    Real costheta = LinAlg3x3::findCosangle(vec1, vec2);
                    Real sintheta = LinAlg3x3::findSinangle(vec1, vec2);

                    factor += costheta/sintheta;
                }

                // find the difference of the edge vector
                Real3 diff;
                diff.fill(0);
                for (int k=0;k<3;k++)
                {
                    diff[k] = vertices[idx1].position_[k] - vertices[idx2].position_[k]; 
                }

                if (idx1 == i)
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
    // obtain the normals of the surfaces 
    mesh_.CalcVertexNormals();

    // Get sum of triangles  --> 
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

    rhs.resize(3, Eigen::VectorXf::Zero(vertexPos_.size()));
    xyz.resize(3, Eigen::VectorXf::Zero(vertexPos_.size()));

    #pragma omp parallel for
    for(int i = 0; i < vertexPos_.size(); ++i)
    {
        // rhs[0][i] = lambdadt_ * Lfactors_[i][0] + vertexPos_[i][0];
        // rhs[1][i] = lambdadt_ * Lfactors_[i][1] + vertexPos_[i][1];
        // rhs[2][i] = lambdadt_ * Lfactors_[i][2] + vertexPos_[i][2];
        rhs[0][i] = vertexPos_[i][0];
        rhs[1][i] = vertexPos_[i][1];
        rhs[2][i] = vertexPos_[i][2];

        xyz[0][i] = vertexPos_[i][0]; 
        xyz[1][i] = vertexPos_[i][1];
        xyz[2][i] = vertexPos_[i][2];
    }

    // solver 
    for (int i=0;i<3;i++)
    {
        xyz[i] = solver_.solveWithGuess(rhs[i],xyz[i]);
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
    // What to do for periodic Meshes?
    /*
        We can make virtual vertices that updates every step.

        For each vertex i, we can check its neighbors, find their distances. If the distance is larger
        than box/2 as usually defined for PBC, we can then make a virtual vertex for vertex i. 
    */
    const auto& vertices = mesh_.getvertices();

    Lfactors_.clear();
    Lfactors_.resize(vertices.size(), {{0,0,0}});

    const auto& neighborIndices = mesh_.getNeighborIndices();
    Real epsilon=1e-5;

    triplets_.clear();
    triplets_buffer_.set_master_object(triplets_);

    // copy the positions of the vertices into a vector
    vertexPos_.clear();
    vertexPos_.resize(vertices.size());
    for (int i=0;i<vertices.size();i++)
    {
        vertexPos_[i] = vertices[i].position_;
    }

    // let's first make all the virtual indices
    std::vector<std::map<int,int>> neighborIndices_corrected;
    // do all the virtual site business if mesh is periodic 
    if (mesh_.isPeriodic())
    {
        auto boxLength = mesh_.getBoxLength();
        for (int i=0;i<vertices.size();i++)
        {
            std::map<int,int> tempMap;
            int m = 0;
            for (auto ind : neighborIndices[i])
            {
                // first check if this is a periodic edge
                Real3 newarr = {};
                bool periodicE = MeshTools::isPeriodicEdge(vertexPos_[ind], vertexPos_[i], newarr , boxLength);
 

                // if it is a periodic Edge, then let's push back on vertexPos vector 
                if (periodicE)
                {
                    tempMap.insert(std::make_pair(ind, vertexPos_.size()));
                    vertexPos_.push_back(newarr);
                    m += 1;
                }
                else
                {
                    tempMap.insert(std::make_pair(ind, ind));
                }
            }
            // push back the temporary index array 
            neighborIndices_corrected.push_back(tempMap);
        }
    }
    else
    {
    // if Mesh is not periodic, then we will just copy the array of neighborIndices 
        // neighborIndices_corrected.insert(neighborIndices_corrected.end(), neighborIndices.begin(), neighborIndices.end());
    }
    std::ofstream ofs;
    ofs.open("test.out");
    std::ofstream ofs2;
    ofs2.open("rhs.out");

    for (int i=0;i<vertexPos_.size();i++)
    {
        ofs2 << vertexPos_[i][0] << " " << vertexPos_[i][1] << " " << vertexPos_[i][2] << "\n";
    }
    ofs2.close();

    ofs << vertexPos_.size() << "\n";


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
                Real3 Lfactor = {};
                std::vector<Real> weights = calculateWeights(i, neighborId, Lfactor);
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
                        auto it = neighborIndices_corrected[i].find(neighborId[j]);
                        ASSERT((it != neighborIndices_corrected[i].end()), "Neighbor " << neighborId[j] << " is not found for vertex " << i);
                        std::cout << "Neighbor " << i << " and " << it -> second << "\n";
                        trip_omp.push_back(triplet(i, it -> second, factor * weights[j]));

                        ofs << i << " " << neighborId[j] << " " << factor * weights[j] << "\n";
                    }

                    for (int j=0;j<3;j++)
                    {
                        Lfactor[j] = Lfactor[j] * factor;
                    }
                    Lfactors_[i] = Lfactor;

                    trip_omp.push_back(triplet(i,i, -factor * w));
                    ofs << i << " " << i << " " << -factor * w << "\n";
                }
            }
        }
    }

    ofs.close();
    ofs2.close();

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

    // set the L matrix 
    L_.setZero();
    L_.resize(vertexPos_.size(), vertexPos_.size());
    L_.setFromTriplets(triplets_.begin(), triplets_.end());

    Eigen::SparseMatrix<Real> I = Eigen::MatrixXf::Identity(vertexPos_.size(), vertexPos_.size()).sparseView();

    L_ = I - lambdadt_*L_;

    // Solve for x, y, z
    // solver_.analyzePattern(L_);
    // solver_.factorize(L_);
    solver_.compute(L_);
}