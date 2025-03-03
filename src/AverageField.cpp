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
    // average the field 
    average();

    // first triangulate the field 
    if (useMC_){
        MarchingCubes_.triangulate_field(field_, *mesh_, isoSurfaceVal_, MCpbc_);
    }
    else{
        SurfaceNets_.triangulate(field_, *mesh_, isoSurfaceVal_);
    }

    // if we are cutting mesh
    if (cut_mesh_){
        MeshTools::CutMesh(*mesh_, cut_vec_, cut_below_);
    }
    std::cout << "Done cut mesh" << "\n";

    // refine the mesh 
    for (int i=0;i<refinementstrat_.size();i++){
        refinementstrat_[i] -> refine(*mesh_);
    }

    // calculate the curvature
    for (int i=0;i<curvatures_.size();i++){
        curvatures_[i]->calculate(*mesh_);
    }
}

void AverageField::printField(std::string name){
    std::ofstream ofs;
    ofs.open(name);
    int Nx,Ny,Nz;
    Nx = field_.getN()[0];
    Ny = field_.getN()[1];
    Nz = field_.getN()[2];
    for (int i=0;i<Nx;i++){
        for (int j=0;j<Ny;j++){
            for (int k=0;k<Nz;k++){
                ofs << i << " " << j << " " << k << " " << field_(i,j,k) << "\n";
            }
        }
    }

    ofs.close();
}


void AverageField::average(){
    auto& fieldVec = field_.accessField();
    Real avgFac = 1.0/simstate_.getTotalFramesToBeCalculated();
    #ifdef CUDA_ENABLED
    // copy data from device 
    _field.memcpyDeviceToHost();
    ASSERT((_field.getTotalSize() == fieldVec.size()), "The GPU and CPU size are different.");
    for (int i=0;i<_field.getTotalSize();i++){
        fieldVec[i] = _field[i] * avgFac;
    }
    #else
    // performing average on the field
    for (int i=0;i<fieldVec.size();i++){
        fieldVec[i] = avgFac*fieldVec[i]; 
    }
    #endif
}