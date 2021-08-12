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
        std::cout << "vertexNormals " << i << " = " << vertexNormals_[i][0] << " " << vertexNormals_[i][1] << " " << vertexNormals_[i][2] << std::endl;
    }

    // for (int i=0;i<vertices_.size();i++)
    // {
    //     vertices_[i].normals_ = vertexNormals_[i];
    // }
}

void AverageCurvature::CalcCurvature()
{
    curvature_.resize(vertices_.size());
    Num_neighbors_.resize(vertices_.size());
    std::fill(Num_neighbors_.begin(), Num_neighbors_.end(),0);
    std::fill(curvature_.begin(), curvature_.end(),1.0);

    // Now we calculate the curvature
    for (int i=0;i<triangles_.size();i++)
    {
        auto& t = triangles_[i];
        for (int j=0;j<3;j++)
        {
            int index1 = t.triangleindices_[j];
            auto& v1 = vertices_[index1];
            for (int k=j+1;k<3;k++)
            {
                int index2 = t.triangleindices_[k];
                auto& v2 = vertices_[index2];

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

                

                Num_neighbors_[index1] += 1;
                Num_neighbors_[index2] += 1;

                curvature_[index1]  *= std::abs(curve_Val);
                curvature_[index2]  *= std::abs(curve_Val);

                if (v2.position_[0] > 6.8)
                {
                    std::cout << "normal2 = " << v2.normals_[0] << " " << v2.normals_[1] << " " << v2.normals_[2] << std::endl;
                    std::cout << "normal1 = " << v1.normals_[0] << " " << v1.normals_[1] << " " << v1.normals_[2] << std::endl;
                    std::cout << "diff = " << diff[0] << " " << diff[1] << " " << diff[2] << std::endl;
                    std::cout << "diffn = " << diffn[0] << " " << diffn[1] << " " << diffn[2] << std::endl;
                    std::cout << "dot_product = " << dot_product << std::endl;
                    std::cout << "sq_diff = " << sq_diff << std::endl;
                    std::cout << "val = " << curve_Val << std::endl;
                    std::cout << "curvature1, index =  " << index1 << " = " << curvature_[index1] << std::endl;
                    std::cout << "curvature2, index = " << index2 << " = " << curvature_[index2] << std::endl; 
                }
            }
        }
    }

    for (int i=0;i<curvature_.size();i++)
    {
        std::cout << "Curvature before = " << curvature_[i] << std::endl;
        std::cout << "Num neighbors = " << Num_neighbors_[i] << std::endl;
        curvature_[i] = std::pow(curvature_[i], 1.0/Num_neighbors_[i]);
        std::cout << "After calc = " << curvature_[i] << std::endl;
    }    
}

void AverageCurvature::finishCalculate()
{
    auto& fieldVec = field_.accessField();
    Real avgFac = 1.0/simstate_.getTotalFramesToBeCalculated();

    // performing average on the field
    for (int i=0;i<fieldVec.size();i++)
    {
        fieldVec[i] = avgFac*fieldVec[i]; 
    }
    
    // calculate the vertices of the system from the averaged density field
    MarchingCubes_.calculate(field_, vertices_, triangles_, isoSurfaceVal_);

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