#include "DensityField.h"

DensityField::DensityField(const DensityFieldInput& input)
:simstate_(input.simstate_)
{
    // Read the dimensions of the system, as in (Nx, Ny, Nz)
    input.pack_.ReadArrayNumber("dimensions", ParameterPack::KeyType::Required, dimensions_);

    // Read the atomgroup names of the system
    input.pack_.ReadString("atomgroup", ParameterPack::KeyType::Required, atomGroupName_);
    addAtomGroup(atomGroupName_);

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

    // Read in the output file
    input.pack_.ReadString("output", ParameterPack::KeyType::Optional, output_name_);

    // calculate the actual cut off
    cutoff_ = n_*sigma_;

    // resize the field 
    field_.resize(dimensions_[0], dimensions_[1], dimensions_[2], x_range_, y_range_, z_range_);

    // calculate the offset Index
    CalcOffsetIndex();

    for (int i=0;i<offsetIndex_.size();i++)
    {
        for (int j=0;j<3;j++)
        {
            std::cout << offsetIndex_[i][j] << " ";
        } 
        std::cout << "\n";
    }

    for (auto it = FieldBuffer_.beginworker();it != FieldBuffer_.endworker();it++)
    {
        it -> resize(dimensions_[0], dimensions_[1], dimensions_[2], x_range_, y_range_, z_range_);
    }
}


void DensityField::CalcOffsetIndex()
{
    // get the differentials in the 3 directions
    Real dx = field_.getdx();
    Real dy = field_.getdy();
    Real dz = field_.getdz();

    int Nx_offset = cutoff_/dx;
    int Ny_offset = cutoff_/dy;
    int Nz_offset = cutoff_/dz;

    for( int i=-Nx_offset/2; i<=Nx_offset/2;i++)
    {
        for (int j=-Ny_offset/2; j<=Ny_offset/2; j++)
        {
            for (int k=-Nz_offset/2; k<=Nz_offset/2;k++)
            {
                index3 id = {{ i,j,k}};
                offsetIndex_.push_back(id);
            }
        }
    }
}

void DensityField::addAtomGroup(std::string& name)
{
    int index = AtomGroups_.size();

    registerAtomGroupID(name, index);

    AtomGroups_.push_back(&simstate_.getAtomGroup(name));
}

void DensityField::registerAtomGroupID(std::string& name, int index)
{
    auto it = MapAtomGroupName2Id_.find(name);

    ASSERT(( it == MapAtomGroupName2Id_.end()), "The name for atomgroup " << name << " is already registered.");

    MapAtomGroupName2Id_.insert(std::make_pair(name, index));
}

int DensityField::getAtomGroupID(std::string& name)
{
    auto it = MapAtomGroupName2Id_.find(name);

    ASSERT(( it != MapAtomGroupName2Id_.end()), "The name for atomgroup " << name << " is already registered.");

    return it -> second;
}

const AtomGroup& DensityField::getAtomGroup(std::string& name)
{
    int ID = getAtomGroupID(name);

    return *AtomGroups_[ID];
}

AtomGroup& DensityField::accessAtomGroup(std::string& name)
{
    int ID = getAtomGroupID(name);

    return const_cast<AtomGroup&>(*AtomGroups_[ID]);
}