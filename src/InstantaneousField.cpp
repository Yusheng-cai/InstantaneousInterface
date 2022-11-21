#include "InstantaneousField.h"

InstantaneousField::InstantaneousField(DensityFieldInput& input)
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

    // triangulate the field --> clears the vertices and triangles of mesh inside this function 
    MarchingCubes_.triangulate_field(field_, *mesh_, isoSurfaceVal_, MCpbc_);

    // once we triangulate the field, we calculate the curvature 
    for (int i=0;i<curvatures_.size();i++){
        curvatures_[i]->calculate(*mesh_);
        const auto& faceC = curvatures_[i]->getFaceCurvature();
        std::vector<Real> faceMeanC(faceC.size(),0.0);

        for (int j=0;j<faceMeanC.size();j++){
            faceMeanC[j] = faceC[j][0];
        }
    }
}