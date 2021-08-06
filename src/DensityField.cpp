#include "DensityField.h"

DensityField::DensityField(const DensityFieldInput& input)
:simstate_(input.simstate_)
{
    // Read the dimensions of the system, as in (Nx, Ny, Nz)
    input.pack_.ReadArrayNumber("dimensions", ParameterPack::KeyType::Required, dimensions_);

    // Read the atomgroup names of the system
    input.pack_.ReadVectorString("atomgroups", ParameterPack::KeyType::Required, atomGroupNames_);

    // Read in the sigma value for the gaussian smoothing function
    input.pack_.ReadNumber("sigma", ParameterPack::KeyType::Required, sigma_);

    // Read in the cut off value (n*sigma_)
    input.pack_.ReadNumber("cutoff", ParameterPack::KeyType::Optional, n_);

    // Read in the bounding box
    input.pack_.ReadString("boundingbox", ParameterPack::KeyType::Required,boundingboxName_);
    bound_box_ = &simstate_.getBoundingBox(boundingboxName_);
    x_range_ = bound_box_->getXrange();
    y_range_ = bound_box_->getYrange();
    z_range_ = bound_box_->getZrange();

    // calculate the actual cut off
    cutoff_ = n_*sigma_;

    // resize the field 
    field_.resize(dimensions_[0], dimensions_[1], dimensions_[2], x_range_, y_range_, z_range_);
}