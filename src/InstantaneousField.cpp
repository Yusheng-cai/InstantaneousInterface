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
}