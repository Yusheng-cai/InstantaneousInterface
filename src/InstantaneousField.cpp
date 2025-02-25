#include "InstantaneousField.h"
namespace DensityFieldRegistry
{
    registry_<InstantaneousField> registerInstantaneousField("InstantaneousField");
}

InstantaneousField::InstantaneousField(const DensityFieldInput& input)
: DensityField(input)
{
    index_=1;
}

void InstantaneousField::calculate()
{
    // zero out the field before each calculation 
    field_.zero();

    // calculate instantaneous field 
    CalculateInstantaneousField();

    #ifdef CUDA_ENABLED
    InstantaneousFieldCopy_cuda();
    #endif

    // triangulate the field --> clears the vertices and triangles of mesh inside this function 
    MarchingCubes_.triangulate_field(field_, *mesh_, isoSurfaceVal_, MCpbc_);

    // if we are cutting mesh
    if (cut_mesh_){
        MeshTools::CutMesh(*mesh_, cut_vec_, cut_below_);
    }

    // once we triangulate the field, we calculate the curvature 
    for (int i=0;i<curvatures_.size();i++){
        curvatures_[i]->calculate(*mesh_);
        const auto& faceC = curvatures_[i]->getFaceCurvature();
        std::vector<Real> faceMeanC(faceC.size(),0.0);

        for (int j=0;j<faceMeanC.size();j++){
            faceMeanC[j] = faceC[j][0];
        }
    }
    std::cout << "Done." << std::endl;
}