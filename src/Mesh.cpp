#include "MarchingCubesWrapper.h"

void Mesh::CalcPerVertexDir()
{
    PerVertexdir1_.resize(vertices_.size());
    PerVertexdir2_.resize(vertices_.size());

    for (int i=0;i<triangles_.size();i++)
    {
        auto& t = triangles_[i];
 
        int index0 = t.triangleindices_[0];
        int index1 = t.triangleindices_[1];
        int index2 = t.triangleindices_[2];

        Real3 diff0;
        Real3 diff1;
        Real3 diff2;

        for (int j=0;j<3;j++)
        {
            diff0[j] = vertices_[index1].position_[j] - vertices_[index0].position_[j];
            diff1[j] = vertices_[index2].position_[j] - vertices_[index1].position_[j];
            diff2[j] = vertices_[index0].position_[j] - vertices_[index2].position_[j];
        }

        PerVertexdir1_[index0] = diff0;
        PerVertexdir1_[index1] = diff1;
        PerVertexdir1_[index2] = diff2;
    }

    for (int i=0;i<PerVertexdir1_.size();i++)
    {
        PerVertexdir1_[i] = LinAlg3x3::CrossProduct(PerVertexdir1_[i], vertices_[i].normals_);

        Real3 B = LinAlg3x3::CrossProduct(vertices_[i].normals_, PerVertexdir2_[i]);
        LinAlg3x3::normalize(B);

        PerVertexdir2_[i] = B;
    } 
}

void Mesh::findVertexNeighbors()
{
    neighborIndices_.clear();
    neighborIndices_.resize(vertices_.size());

    NumNeighbors_.clear();
    NumNeighbors_.resize(vertices_.size());

    for (int i=0;i<triangles_.size();i++)
    {
        auto& t = triangles_[i];

        for (int j=0;j<3;j++)
        {
            int index1 = t.triangleindices_[j];
            for (int k=0;k<3;k++)
            {
                if (j != k)
                {
                    int index2 = t.triangleindices_[k];

                    auto it = std::find(neighborIndices_[index1].begin(), neighborIndices_[index1].end(), index2);

                    if (it == neighborIndices_[index1].end())
                    {
                        neighborIndices_[index1].push_back(index2);
                        NumNeighbors_[index1] += 1;
                    }
                }
            }
        }
    }
}

void Mesh::findTriangleIndices()
{
    VertexTriangleIndices_.clear();
    VertexTriangleIndices_.resize(vertices_.size());

    for (int i=0;i<triangles_.size();i++)
    {
        auto& t = triangles_[i];

        for (int j=0; j<3;j++)
        {
            int index = t.triangleindices_[j];

            VertexTriangleIndices_[index].push_back(i);
        }
    }
}

void Mesh::CalcTriangleAreaAndFacetNormals()
{
    triangleArea_.clear();
    triangleArea_.resize(triangles_.size());

    facetNormals_.clear();
    facetNormals_.resize(triangles_.size());

    triangleArea_.clear();
    triangleArea_.resize(triangles_.size());

    // calculate the area of the triangles as well as the normals of the faces of triangles
    for (int i=0;i<triangles_.size();i++)
    {
        auto& t = triangles_[i];
        int index1 = t.triangleindices_[0];
        int index2 = t.triangleindices_[1];
        int index3 = t.triangleindices_[2];

        Real3 diff1;
        Real3 diff2;
        diff1.fill(0);
        diff2.fill(0);

        for (int i=0;i<3;i++)
        {
            diff1[i] = vertices_[index2].position_[i] - vertices_[index1].position_[i];
            diff2[i] = vertices_[index1].position_[i] - vertices_[index3].position_[i];
        }

        Real3 crossProduct = LinAlg3x3::CrossProduct(diff1, diff2);
        Real norm = LinAlg3x3::norm(crossProduct);
        Real Area = norm*0.5;
        Real3 normal_;
        for (int i=0;i<3;i++)
        {
            normal_[i] = crossProduct[i]/norm;
        }

        facetNormals_[i] = normal_;
        triangleArea_[i] = Area;
    }
}

void Mesh::CalcVertexNormals()
{
    // calculate vertex normals
    vertexNormals_.resize(vertices_.size());
    Real3 zeroArr = {{0,0,0}};
    std::fill(vertexNormals_.begin(), vertexNormals_.end(), zeroArr);

    // calculate the normals of each vertices
    for (int i=0;i<triangles_.size();i++)
    {
        auto& t = triangles_[i];
        Real area = triangleArea_[i];
        Real3 normals = facetNormals_[i];

        for (int j=0;j<3;j++)
        {
            int index = t.triangleindices_[j];
            auto& vNorm = vertexNormals_[index];

            for (int k=0;k<3;k++)
            {
                vNorm[k] += area*normals[k];
            }
        }
    }


    for (int i=0;i<vertexNormals_.size();i++)
    {
        Real norm = LinAlg3x3::norm(vertexNormals_[i]);

        for (int j=0;j<3;j++)
        {
            vertexNormals_[i][j] = vertexNormals_[i][j]/norm;
        }
    }
}