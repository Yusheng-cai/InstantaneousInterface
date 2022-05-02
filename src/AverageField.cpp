#include "AverageField.h"

namespace DensityFieldRegistry
{
    registry_<AverageField> registerAverageField("averagefield");
}

AverageField::AverageField(const DensityFieldInput& input)
:DensityField(input)
{
    outputs_.registerOutputFunc("field", [this](std::string name) -> void { this -> printField(name);});
}

void AverageField::calculate()
{
    CalculateInstantaneousField();
}

void AverageField::finishCalculate()
{
    auto& fieldVec = field_.accessField();

    Real avgFac = 1.0/simstate_.getTotalFramesToBeCalculated();

    // performing average on the field
    for (int i=0;i<fieldVec.size();i++)
    {
        fieldVec[i] = avgFac*fieldVec[i]; 
    }

    // first triangulate the field 
    MarchingCubes_.triangulate_field(field_, *mesh_, isoSurfaceVal_, MCpbc_);

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