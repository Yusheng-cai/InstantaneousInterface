#include "MeshLoopRefine.h"

namespace MeshRefineStrategyFactory
{
    registry_<MeshLoopRefine> registerLoopRefinement("loop");
}

MeshLoopRefine::MeshLoopRefine(MeshRefineStrategyInput& input)
:MeshRefineStrategy(input)
{
}

void MeshLoopRefine::refine()
{
    // calculate the even points of the loop refinement 
    calculateEvenPoints();

    // calculate the odd points of the loop refinement 
    calculateOddPoints();

    // now calculate the triangles
    calculateTriangles();

    // now update the mesh 
    updateMesh();
}

void MeshLoopRefine::updateMesh()
{
    auto& vertices = mesh_.accessvertices();
    auto& triangles = mesh_.accesstriangles();
    vertices.clear();
    triangles.clear();

    // insert vertices
    vertices.insert(vertices.end(),Points_.begin(), Points_.end());
    triangles.insert(triangles.end(), triangles_.begin(), triangles_.end());

    // calculate triangles faces
    mesh_.CalcTriangleAreaAndFacetNormals();
}

void MeshLoopRefine::calculateTriangles()
{
    // each face is split into 4 smaller triangles 
    const auto& oldTriangles = mesh_.gettriangles();

    for (int i=0;i<oldTriangles.size();i++)
    {
        // the edges of a triangle
        auto edges = oldTriangles[i].edges_;

        // the indices of the points in the old face 
        for (int j=0;j<3;j++)
        {
            std::vector<int> indices_;
            int id = oldTriangles[i].triangleindices_[j];
            indices_.push_back(id);

            for (auto e:edges)
            {
                auto it = MapEdgesToVertex_.find(e);
                ASSERT((it != MapEdgesToVertex_.end()), "The edge is not found.");

                if (id == e.vertex1_.index)
                {
                    indices_.push_back(it -> second.index);
                }

                if (id == e.vertex2_.index)
                {
                    indices_.push_back(it -> second.index);
                }
            }

            ASSERT((indices_.size() == 3), "The triangular indices cannot exceed 3 while it is " << indices_.size());
            index3 idd_;

            for (int k=0;k<3;k++)
            {
                idd_[k] = indices_[k];
            }

            triangle_indices_.push_back(idd_);
        }    

        index3 id_;
        // Now we connect all the inner triangles
        for (int j=0;j<3;j++)
        {
            edge e = edges[j];

            auto it = MapEdgesToVertex_.find(e);

            ASSERT((it != MapEdgesToVertex_.end()), "The edge with e.vertex1 = " << e.vertex1_.index << " and e.vertex2 = " << e.vertex2_.index << \
            " is not in Map from edges to vertex.");

            id_[j] = it -> second.index;
        }

        triangle_indices_.push_back(id_);
    }

    // Now we calculate the normal vector of the triangles and construct the triangles 
    std::vector<Real> TotalAreaVertices_(Points_.size(),0.0);
    triangles_.resize(triangle_indices_.size());
    for (int i=0;i<triangles_.size();i++)
    {
        Real3 edge1;
        Real3 edge2;
        auto& t = triangle_indices_[i];

        for (int j=0;j<3;j++)
        {
            edge1[j] = Points_[t[0]].position_[j] - Points_[t[1]].position_[j];
            edge2[j] = Points_[t[1]].position_[j] - Points_[t[2]].position_[j];
        }

        Real3 normal = LinAlg3x3::CrossProduct(edge1, edge2);
        Real area = LinAlg3x3::norm(normal);
        LinAlg3x3::normalize(normal);

        for (int j=0;j<3;j++)
        {
            for (int k=0;k<3;k++)
            {
                Points_[t[j]].normals_[k] += area * normal[k];
            }

            TotalAreaVertices_[t[j]] += area;
        }
    }

    for (int i=0;i<Points_.size();i++)
    {
        Real TotalArea = TotalAreaVertices_[i];

        for (int j=0;j<3;j++)
        {
            Points_[i].normals_[j] /= TotalArea;
        }
    }

    for (int i=0;i<triangles_.size();i++)
    {
        for (int j=0;j<3;j++)
        {
            int vertexIndex = triangle_indices_[i][j];
            triangles_[i].vertices_[j] = Points_[vertexIndex];
            triangles_[i].triangleindices_[j] = vertexIndex;
        }
        triangles_[i].edges_[0].vertex1_ = triangles_[i].vertices_[0];
        triangles_[i].edges_[0].vertex2_ = triangles_[i].vertices_[1];

        triangles_[i].edges_[1].vertex1_ = triangles_[i].vertices_[1];
        triangles_[i].edges_[1].vertex2_ = triangles_[i].vertices_[2];

        triangles_[i].edges_[2].vertex1_ = triangles_[i].vertices_[2];
        triangles_[i].edges_[2].vertex2_ = triangles_[i].vertices_[0];
    }
}

void MeshLoopRefine::calculateOddPoints()
{
    // calculate the map that maps edges to Faces
    mesh_.MapEdgeToFaces();
    const auto& map = mesh_.getMapEdgeToFace();
    const auto& triangles = mesh_.gettriangles();
    const auto& vertices = mesh_.getvertices();
    
    // first calculate the odd points
    for (int i=0;i<map.size(); i++)
    {
        auto it = map.begin();
        std::advance(it,i);

        // find the edge
        edge e = it -> first;

        int NumFaces = it -> second.size();
        auto& faceIndex = it -> second;

        ASSERT((NumFaces == 1 || NumFaces == 2), "Number of faces of shared by a single edge needs to be either 1 or 2.");

        // first take care of the boundary cases
        if (NumFaces == 1)
        {
            Real3 pos;
            pos.fill(0);

            vertex v;
            for (int i=0;i<3;i++)
            {
                pos[i] = 0.5*(e.vertex1_.position_[i] + e.vertex2_.position_[i]);
            }
            v.position_ = pos;
            v.index     = Points_.size();
            v.normals_  = {{0,0,0}};
            Points_.push_back(v);

            // map edges to their respective vertex 
            auto it = MapEdgesToVertex_.find(e);
            ASSERT((it == MapEdgesToVertex_.end()), "The edge is already registered, each edge should only have 1 point associated with it.");
            MapEdgesToVertex_.insert(std::make_pair(e, v));
        }
        else
        {
            vertex v;
            Real3 pos;
            pos.fill(0);

            ASSERT((faceIndex.size() == 2), "Each middle edges can only be shared by 2 faces, however calculation suggests that it is shared by " \
             << faceIndex.size());

            // the triangle pair 
            auto t1 = triangles[faceIndex[0]];
            auto t2 = triangles[faceIndex[1]];
            std::array<triangle,2> trPair = {{t1,t2}};

            std::array<int,2> thisEdge = {{e.vertex1_.index, e.vertex2_.index}};
            std::vector<int> OtherEdge;

            #ifdef MY_DEBUG
            std::cout << "New edge" << std::endl;
            for (int i=0;i<2;i++)
            {
                auto tIndices = trPair[i].triangleindices_;
                for (int j=0;j<3;j++)
                {
                    std::cout << tIndices[j] << "\t";
                }
            }
            std::cout << "\n";
            #endif 

            for (int i=0;i<2;i++)
            {
                auto tIndices = trPair[i].triangleindices_;
                for (int j=0;j<3;j++)
                {
                    // first check if the indices is already in other edge
                    bool NotinOtherEdge = (std::find(OtherEdge.begin(), OtherEdge.end(), tIndices[j])) == OtherEdge.end();

                    // then check if the indices is in this Edge
                    bool NotinThisEdge = (std::find(thisEdge.begin(), thisEdge.end(), tIndices[j])) == thisEdge.end();

                    if (NotinThisEdge && NotinOtherEdge)
                    {
                        OtherEdge.push_back(tIndices[j]);
                    }
                }
            }

            // we assert the other edge has only 2 indices
            ASSERT((OtherEdge.size() == 2), "The other edge must have 2 vertices while it has " << OtherEdge.size());

            // then we apply the formula for the odd points
            for (int i=0;i<3;i++)
            {
                Real thisEdgeval = vertices[thisEdge[0]].position_[i] + vertices[thisEdge[1]].position_[i];
                Real otherEdgeval = vertices[OtherEdge[0]].position_[i] + vertices[OtherEdge[1]].position_[i];

                pos[i] += 0.125*otherEdgeval + 0.375*thisEdgeval;
            }

            v.position_ = pos;
            v.index = Points_.size();
            v.normals_ = {{0,0,0}};
            Points_.push_back(v);

            // Map edges to vertex 
            auto it2 = MapEdgesToVertex_.find(e);
            ASSERT((it2 == MapEdgesToVertex_.end()), "The edge is already registered, each edge should only have 1 point associated with it.");
            MapEdgesToVertex_.insert(std::make_pair(e, v));
        }
    }
}

void MeshLoopRefine::calculateEvenPoints()
{
    // first we find the neighbors of the vertex
    mesh_.findVertexNeighbors();

    // find the boundary points
    mesh_.findBoundaryVertices();

    const auto& neighbors = mesh_.getNeighborIndices();
    const auto& vertices  = mesh_.getvertices();
    const auto& bVertTobEdge = mesh_.getMapBVertexToBEdges();

    // resize the evenpoints
    int numPoints = neighbors.size();

    for (int i=0;i<numPoints;i++)
    {
        int lenNeighbors = neighbors[i].size();
        Real beta = this->calculateBeta(lenNeighbors);

        Real3 posI = vertices[i].position_;
        Real3 sum_ = {};

        vertex v;

        // see if this point is a boundary point
        auto it = bVertTobEdge.find(vertices[i]);

        // if the point is not a boundary vertices 
        if (it == bVertTobEdge.end())
        {
            for (int j=0;j<lenNeighbors;j++)
            {
                Real3 posN = vertices[neighbors[i][j]].position_;
                for (int k=0;k<3;k++)
                {
                    sum_[k] += posN[k];
                }
            }

            for (int j=0;j<3;j++)
            {
                sum_[j] = (1 - lenNeighbors*beta)*posI[j] + beta * sum_[j];
            }
        }
        else
        {
            // if the point is a boundary vertex 
            std::vector<edge> edges = it -> second;

            // assert that there must be 2 edges
            ASSERT((edges.size() == 2), "The number of edges for a boundary point must be 2 while it is " << edges.size());

            std::vector<Real3> OtherPoints;

            for (auto e:edges)
            {
                if (e.vertex1_.index != i)
                {
                    OtherPoints.push_back(e.vertex1_.position_);
                }

                if (e.vertex2_.index != i)
                {
                    OtherPoints.push_back(e.vertex2_.position_);
                }
            }

            ASSERT((OtherPoints.size() == 2), "It must be 2.");

            for (int j=0;j<3;j++)
            {
                sum_[j] += 1.0/8.0*(OtherPoints[0][j] + OtherPoints[1][j]) + 3.0/4.0*posI[j];
            }
        }

        v.position_ = sum_;
        v.index = Points_.size();
        v.normals_ = {{0,0,0}};

        // we will calculate normals later 

        // push back on the even points
        Points_.push_back(v);
    }
}

MeshLoopRefine::Real MeshLoopRefine::calculateBeta(int n)
{
    if (n == 3)
    {
        return 3.0/16.0;
    }
    else
    {
        Real factor = 1.0/(Real)n;
        Real l = 3.0/8.0 + 0.25*std::cos(6.28318531/(Real)n);
        Real l2 = l*l;

        Real beta = factor * (5.0/8.0 - 0.25*l2);

        return beta;
    }
}
