#include "AverageCurvature.h"

namespace DensityFieldRegistry
{
    registry_<AverageCurvature> registerAverageCurvature("averagecurvature");
}

AverageCurvature::AverageCurvature(const DensityFieldInput& input)
:AverageField(input)
{
    bool readCurve = input.pack_.ReadString("curvature_output", ParameterPack::KeyType::Optional, curvatureOutputName_);

    if (readCurve)
    {
        curvatureofs_.open(curvatureOutputName_);
        ASSERT((curvatureofs_.is_open()), "The file with name " << curvatureOutputName_ << " is not opened.");
    }
}

void AverageCurvature::CalcNeighborsOfVertex()
{
    Num_neighbors_.clear();
    neighbor_indices_.clear();
    vertex_triangle_indices_.clear();

    Num_neighbors_.resize(vertices_.size());
    vertex_triangle_indices_.resize(vertices_.size());
    std::fill(Num_neighbors_.begin(), Num_neighbors_.end(), 0);
    neighbor_indices_.resize(vertices_.size());

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

                    auto it = std::find(neighbor_indices_[index1].begin(), neighbor_indices_[index1].end(), index2);

                    if (it == neighbor_indices_[index1].end())
                    {
                        neighbor_indices_[index1].push_back(index2);
                        Num_neighbors_[index1] += 1;

                        vertex_triangle_indices_[index1].push_back(i);
                    }
                }
            }
        }
    }
}

void AverageCurvature::CalcTriangleAreaAndFacetNormals()
{
    TriangleAreas_.resize(triangles_.size());
    TriangleNormals_.resize(triangles_.size());

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

        TriangleNormals_[i] = normal_;

        TriangleAreas_[i] = Area;
    }
}

void AverageCurvature::CalcVertexNormals()
{
    // calculate vertex normals
    vertexNormals_.resize(vertices_.size());
    Real3 zeroArr = {{0,0,0}};
    std::fill(vertexNormals_.begin(), vertexNormals_.end(), zeroArr);

    // calculate the normals of each vertices
    for (int i=0;i<triangles_.size();i++)
    {
        auto& t = triangles_[i];
        Real area = TriangleAreas_[i];
        Real3 normals = TriangleNormals_[i];

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

void AverageCurvature::CalcCurvature()
{
    curvature_.clear();
    curvature_.resize(vertices_.size());
    std::fill(curvature_.begin(), curvature_.end(),0.0);


    neighborAreaTotlist_.clear();
    neighborAreaTotlist_.resize(vertices_.size());
    std::fill(neighborAreaTotlist_.begin(), neighborAreaTotlist_.end(),0.0);

    // Now we calculate the curvature
    for (int i=0;i<vertices_.size();i++)
    {
        auto& v1 = vertices_[i];
        auto& neighbors = neighbor_indices_[i];
        auto& triangleInd = vertex_triangle_indices_[i];

        ASSERT(( triangleInd.size() == neighbors.size()), "The neighbor triangle indices does not match the number of neighbors.");

        Real neighborAreaTot = 0.0;
        for (int j=0;j<neighbors.size();j++)
        {
            int id = triangleInd[j];
            neighborAreaTot += TriangleAreas_[id];
        }
        neighborAreaTotlist_[i] = neighborAreaTot;

        for (int j=0;j<neighbors.size();j++)
        { 
            auto& v2 = vertices_[neighbors[j]];
            auto& triangleI = triangleInd[j];

            Real3 diff;
            Real3 diffn;
            Real sq_diff=0.0;
            Real dot_product = 0.0;
            Real curve_Val;
            diff.fill(0);
            diffn.fill(0);
            
            for(int m=0;m<3;m++)
            {
                diff[m] = v2.position_[m] - v1.position_[m];
                diffn[m] = v2.normals_[m] - v1.normals_[m];
            }

            Real diffnnorm = LinAlg3x3::norm(diffn); 

            for (int m=0;m<3;m++)
            {
                sq_diff += diff[m]*diff[m];
            }

            for (int m=0;m<3;m++)
            {
                dot_product += diff[m]*diffn[m];
            }

            curve_Val = dot_product/sq_diff;

            std::cout << "For vertex " << i << " triangle area " << j << "="<< TriangleAreas_[triangleI];
            std::cout << "Total area = " << neighborAreaTotlist_[i] << std::endl;
            curvature_[i]  += std::abs(curve_Val) * TriangleAreas_[triangleI]/neighborAreaTotlist_[i];
        }
    }

    for (int i=0;i<curvature_.size();i++)
    { 
        curvature_[i] = std::exp(curvature_[i]);
    }    
}

void AverageCurvature::finishCalculate()
{
    auto& fieldVec = field_.accessField();
    std::cout << "The number of frames to be averaged = " << simstate_.getTotalFramesToBeCalculated() << std::endl;
    Real avgFac = 1.0/simstate_.getTotalFramesToBeCalculated();

    // performing average on the field
    for (int i=0;i<fieldVec.size();i++)
    {
        fieldVec[i] = avgFac*fieldVec[i]; 
    }
    
    // calculate the vertices of the system from the averaged density field
    MarchingCubes_.calculate(field_, vertices_, triangles_, isoSurfaceVal_);

    CalcNeighborsOfVertex();
    CalcTriangleAreaAndFacetNormals();
    CalcVertexNormals();
    CalcCurvature(); 
}

void AverageCurvature::printCurvature()
{
    if (curvatureofs_.is_open())
    {
        for (int i=0;i<curvature_.size();i++)
        {
            curvatureofs_ << curvature_[i];
            curvatureofs_ << "\n"; 
        }
        curvatureofs_.close();
    }
}

void AverageCurvature::printFinalOutput()
{
    printField();
    printVertices();
    printTriangleIndices();
    printNormals();
    printCurvature();
}