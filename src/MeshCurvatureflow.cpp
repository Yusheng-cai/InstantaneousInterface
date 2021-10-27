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
}

void MeshCurvatureflow::refine()
{
    mesh_.MapEdgeToFaces();

    // find the boundary vertices
    mesh_.findBoundaryVertices();

    // find map from vertex to the face indices 
    mesh_.MapVertexToFaces();

    // find number of vertices 
    const auto& v = mesh_.getvertices();
    newVertices_.insert(newVertices_.end(), v.begin(), v.end());
    TotalArea_.resize(v.size());


    for (int i=0;i<numIterations_;i++)
    {
        refineStep();
    }
}

void MeshCurvatureflow::refineStep()
{
    std::cout << "refinestep" << std::endl;
    mesh_.CalcTriangleAreaAndFacetNormals();

    const auto& vertexIndicesToface = mesh_.getMapVertexToFace();
    const auto& triangleArea = mesh_.getTriangleArea();
    const auto& vertices = mesh_.getvertices();
    const auto& triangles= mesh_.gettriangles();

    std::fill(TotalArea_.begin(), TotalArea_.end(), 0.0);
    std::cout << "Totalarea.size = " << TotalArea_.size() << std::endl;
    std::cout << "vertexindicestoface.size = " << vertexIndicesToface.size() << std::endl;

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

                    for (int m=0;m<3;m++)
                    {
                        vec1[m] = vertices[index1].position_[m] - vertices[pidx].position_[m];
                        vec2[m] = vertices[pidx].position_[m] - vertices[index2].position_[m];
                    }

                    Real costheta = LinAlg3x3::findCosangle(vec1, vec2);
                    Real sintheta = LinAlg3x3::findSinangle(vec1, vec2);

                    std::cout << "Sum of cos2theta and sin2theta is " << costheta*costheta + sintheta*sintheta << std::endl;
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
                newVertices_[i].position_[j] = vertices[i].position_[j] + 1.0/(factor_sum) * lambdadt_ * step[j];
            }
        }
    }

    auto& vert = mesh_.accessvertices();
    vert.clear();
    vert.insert(vert.end(), newVertices_.begin(), newVertices_.end());
    mesh_.update();
    std::cout << "Mesh updated." << std::endl;
}