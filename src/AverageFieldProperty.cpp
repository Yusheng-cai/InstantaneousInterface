#include "AverageFieldProperty.h"

namespace DensityFieldRegistry
{
    registry_<AverageFieldProperty> registerAverageFieldProerty("averagefieldproperty");
}

AverageFieldProperty::AverageFieldProperty(const DensityFieldInput& input)
:DensityField(input){
    outputs_.registerOutputFunc("field", [this](std::string name) -> void { this -> printField(name);});
    pack_.ReadString("property_file", ParameterPack::KeyType::Required, property_filename_);

    // read the properties
    StringTools::ReadTabulatedData(property_filename_, properties_);
    ASSERT((properties_.size() >= simstate_.getTotalFrames()), "Properties size is less than total frames to be calculated.");
}

void AverageFieldProperty::calculate(){
    int frame_num = simstate_.getFrame();
    ASSERT((frame_num < properties_.size()), "Out of bounds for properties.");

    CalculateInstantaneousFieldProperty(properties_[frame_num]);
}

void AverageFieldProperty::finishCalculate()
{
    // average the field 
    average();

    // first triangulate the field 
    if (useMC_){
        MarchingCubes_.triangulate_field(property_field_, *mesh_, isoSurfaceVal_, MCpbc_);
    }
    else{
        SurfaceNets_.triangulate(property_field_, *mesh_, isoSurfaceVal_);
    }

    // if we are cutting mesh
    if (cut_mesh_){
        MeshTools::CutMesh(*mesh_, cut_vec_);
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

void AverageFieldProperty::printField(std::string name){
    std::ofstream ofs;
    ofs.open(name);
    for (int i=0;i<property_field_.accessField().size();i++) {
        ofs << property_field_.accessField()[i];
        ofs << " ";
    }

    ofs.close();
}


void AverageFieldProperty::average(){
    auto& fieldVec = field_.accessField();
    auto& propertyVec = property_field_.accessField();

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
    for (int i=0;i<propertyVec.size();i++){
        propertyVec[i] = propertyVec[i] * avgFac;
    }
    #endif
}

