#include "AverageField.h"

namespace DensityFieldRegistry
{
    registry_<AverageField> registerAverageField("averagefield");
}

AverageField::AverageField(const DensityFieldInput& input)
:DensityField(input)
{
    outputs_.registerOutputFunc("field", [this](std::string name) -> void { this -> printField(name);});
    outputs_.registerOutputFunc("triangle", [this](std::string name) -> void {this -> printTriangleIndices(name);});
    outputs_.registerOutputFunc("vertex", [this](std::string name) -> void {this -> printVertices(name);});
    outputs_.registerOutputFunc("normal", [this](std::string name) -> void { this -> printNormals(name);});
}

void AverageField::calculate()
{
    CalculateInstantaneousInterface();
}

void AverageField::finishCalculate()
{
    auto& fieldVec = field_.accessField();

    #ifdef MY_DEBUG
    std::cout << "Total frames to be calculated = " << simstate_.getTotalFramesToBeCalculated() << std::endl;
    #endif 

    Real avgFac = 1.0/simstate_.getTotalFramesToBeCalculated();

    // performing average on the field
    for (int i=0;i<fieldVec.size();i++)
    {
        fieldVec[i] = avgFac*fieldVec[i]; 
    }

    MarchingCubes_.triangulate_field(field_, *mesh_, isoSurfaceVal_, pbc_);

    // refine the mesh 
    for (int i=0;i<refinementstrat_.size();i++)
    {
        refinementstrat_[i] -> refine(*mesh_);
    }

    // calculate the curvature
    for (int i=0;i<curvatures_.size();i++)
    {
        curvatures_[i]->calculate(*mesh_);
    }
}

void AverageField::printOutputIfOnStep()
{}

void AverageField::printField(std::string name)
{
    std::ofstream ofs_;
    ofs_.open(name);
    for (int i=0;i<field_.accessField().size();i++) 
    {
        ofs_ << field_.accessField()[i];

        ofs_ << " ";
    }

    ofs_.close();
}

void AverageField::printVertices(std::string name)
{
    std::ofstream ofs_;
    ofs_.open(name);

    const auto& vertices_ = mesh_ -> getvertices();
    for (int i=0;i<vertices_.size();i++)
    {
        for (int j=0;j<3;j++)
        {
            ofs_ << vertices_[i].position_[j];

            ofs_ << " ";
        }

        ofs_ << "\n";
    }

    ofs_.close();
}

void AverageField::printTriangleIndices(std::string name)
{
    std::ofstream ofs_;
    ofs_.open(name);

    const auto& triangles_ = mesh_ -> gettriangles();
    for (int i=0;i<triangles_.size();i++)
    {
        for (int j=0;j<3;j++)
        {
            ofs_ << triangles_[i].triangleindices_[j];

            ofs_ << " ";
        }
        ofs_ << "\n";
    }

    ofs_.close();
}

void AverageField::printNormals(std::string name)
{
    std::ofstream ofs_;
    ofs_.open(name);

    const auto& vertices_ = mesh_ -> getvertices();
    for (int i=0;i<vertices_.size();i++)
    {
        for (int j=0;j<3;j++)
        {
            ofs_ << vertices_[i].normals_[j];

            ofs_ << " ";
        }

        ofs_ << "\n";
    }
    ofs_.close();
}