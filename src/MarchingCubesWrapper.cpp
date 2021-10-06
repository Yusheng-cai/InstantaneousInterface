#include "MarchingCubesWrapper.h"

void MarchingCubesWrapper::calculate(Field& field_, Mesh& mesh_,Real isoSurfaceVal)
{
    index3 index = field_.getN();
    Real3 spacing= field_.getSpacing();
    std::vector<float> vertex_data;
    std::vector<int> triangle_indices;
    int NumberVertex = 0;
    int NumberTriangle = 0;

    marchingCubesalgo_.Triangulate(field_.data(), vertex_data, triangle_indices, index.data(), spacing.data(), &NumberVertex,\
        &NumberTriangle, isoSurfaceVal, 0);
    
    Real3 Length = field_.getLength();

    auto& vertices = mesh_.accessvertices();
    auto& triangles= mesh_.accesstriangles();

    vertices.clear();
    triangles.clear();

    vertices.resize(NumberVertex);
    triangles.resize(NumberTriangle);


    for (int i=0;i<vertices.size();i++)
    {
        vertex vert;
        Real normal_sq = 0.0;
        for (int j=0;j<3;j++)
        {
            // construct the vertices
            vert.position_[j] = vertex_data[6*i+j] + 0.5*Length[j];
            vert.normals_[j]  = vertex_data[6*i+3+j];
            normal_sq += vert.normals_[j]*vert.normals_[j];
        }

        // each vertex in the triangulated mesh has an index
        vert.index = i;

        Real norm = std::sqrt(normal_sq);

        for (int j=0;j<3;j++)
        {
            vert.normals_[j] = vert.normals_[j]/norm;
        }

        vertices[i] = vert;
    }

    for (int i=0;i<triangles.size();i++)
    {
        triangle t; 

        for (int j=0;j<3;j++)
        {
            t.triangleindices_[j] = triangle_indices[3*i + j]; 
            int index = t.triangleindices_[j];

            t.vertices_[j] = vertices[index];
        }

        // construct the edges of each triangle
        t.edges_[0].vertex1_ = t.vertices_[0];
        t.edges_[0].vertex2_ = t.vertices_[1];

        t.edges_[1].vertex1_ = t.vertices_[1];
        t.edges_[1].vertex2_ = t.vertices_[2];

        t.edges_[2].vertex1_ = t.vertices_[2];
        t.edges_[2].vertex2_ = t.vertices_[0];

        triangles[i] = t;
    }
    
    // calculate the triangular area and facet normals
    mesh_.CalcTriangleAreaAndFacetNormals();
    // mesh_.CalcVertexNormals();
}