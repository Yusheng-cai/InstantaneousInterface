#include "MarchingCubesWrapper.h"

void MarchingCubesWrapper::calculate(Field& field_, std::vector<vertex>& vertices, std::vector<triangle>& triangles,Real isoSurfaceVal)
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
            vert.position_[j] = vertex_data[6*i+j] + 0.5*Length[j];
            vert.normals_[j]  = vertex_data[6*i+3+j];
            normal_sq += vert.normals_[j]*vert.normals_[j];
        }

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
        }

        triangles[i] = t;
    }
}