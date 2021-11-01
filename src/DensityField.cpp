#include "DensityField.h"

DensityField::DensityField(const DensityFieldInput& input)
:simstate_(input.simstate_), pack_(input.pack_)
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

    // Read in the isoSurface value
    input.pack_.ReadNumber("isosurfacevalue", ParameterPack::KeyType::Required, isoSurfaceVal_);

    // Read in the bounding box
    input.pack_.ReadString("boundingbox", ParameterPack::KeyType::Required,boundingboxName_);
    bound_box_ = &simstate_.getBoundingBox(boundingboxName_);
    x_range_ = bound_box_->getXrange();
    y_range_ = bound_box_->getYrange();
    z_range_ = bound_box_->getZrange();

    // read in the output names 
    input.pack_.ReadVectorString("outputs", ParameterPack::KeyType::Optional, OutputNames_);

    // read in the output file names 
    input.pack_.ReadVectorString("outputFiles", ParameterPack::KeyType::Optional, OutputFileNames_);
    ASSERT(( OutputFileNames_.size() == OutputNames_.size()), "The number of outputs does not agree with number of \
    outputs files.");

    // calculate the actual cut off
    cutoff_ = n_*sigma_;

    // resize the field 
    field_.resize(dimensions_[0], dimensions_[1], dimensions_[2], x_range_, y_range_, z_range_);

    // calculate the offset Index
    CalcOffsetIndex();

    for (auto it = FieldBuffer_.beginworker();it != FieldBuffer_.endworker();it++)
    {
        it -> resize(dimensions_[0], dimensions_[1], dimensions_[2], x_range_, y_range_, z_range_);
    }

    initializeMesh();
    initializeCurvature();
}

void DensityField::initializeCurvature()
{
    auto pack = pack_.findParamPacks("curvature", ParameterPack::KeyType::Optional);

    if( pack.size() != 0)
    {
        for (int i=0;i<pack.size();i++)
        {
            std::string curvatureType;
            auto& p = pack[i];

            p -> ReadString("type", ParameterPack::KeyType::Required, curvatureType);

            CurvatureInput input = { const_cast<ParameterPack&>(*(p)), *mesh_ };
            CurvaturePtr ptr = CurvaturePtr(CurvatureRegistry::Factory::instance().create(curvatureType, input));

            curvatures_.push_back(std::move(ptr)); 
        }
    }
}

void DensityField::printFinalOutput()
{
    for (int i=0;i<OutputNames_.size();i++)
    {
        outputs_.getOutputFuncByName(OutputNames_[i])(OutputFileNames_[i]);
    }

    for (int i=0;i<curvatures_.size();i++)
    {
        curvatures_[i]->printOutput();
    }

    mesh_ -> print();
}

void DensityField::CalcOffsetIndex()
{
    offsetIndex_.clear();

    // get the differentials in the 3 directions
    Real dx = field_.getdx();
    Real dy = field_.getdy();
    Real dz = field_.getdz();

    int Nx_offset = std::round(cutoff_/dx);
    int Ny_offset = std::round(cutoff_/dy);
    int Nz_offset = std::round(cutoff_/dz);

    for( int i=-Nx_offset; i<=Nx_offset;i++)
    {
        for (int j=-Ny_offset; j<=Ny_offset; j++)
        {
            for (int k=-Nz_offset; k<=Nz_offset;k++)
            {
                index3 id = {{ i,j,k}};
                offsetIndex_.push_back(id);
            }
        }
    }
    #ifdef MY_DEBUG
    std::cout << "Printing offset indices with cutoff = " << cutoff_ << " Nxoffset = " << Nx_offset << " , Nyoffset = " << Ny_offset << " Nzoffset = " << Nz_offset << std::endl;
    for (int i=0;i<offsetIndex_.size();i++)
    {
        std::cout << offsetIndex_[i][0] << " " << offsetIndex_[i][1] << " " << offsetIndex_[i][2] << std::endl;
    }
    #endif 
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

void DensityField::initializeMesh()
{
    auto meshP = pack_.findParamPack("Mesh", ParameterPack::KeyType::Optional);
    mesh_ = Meshptr(new Mesh(meshP));
}

void DensityField::findAtomsIndicesInBoundingBox()
{
    auto atomgroup = getAtomGroup(atomGroupName_);
    auto& atoms = atomgroup.getAtoms();

    AtomIndicesInside_.clear();
    AtomIndicesBuffer_.clearBuffer();
    AtomIndicesBuffer_.set_master_object(AtomIndicesInside_);

    #pragma omp parallel
    {
        auto& indices_buffer = AtomIndicesBuffer_.access_buffer_by_id();
        #pragma omp for
        for(int i=0;i<atoms.size();i++)
        {
            if (bound_box_ -> isInside(atoms[i].position))
            {
                indices_buffer.push_back(i);
            }
        }
    }

    int size = AtomIndicesInside_.size();
    for (auto it = AtomIndicesBuffer_.beginworker();it != AtomIndicesBuffer_.endworker();it++)
    {
        size += it -> size();
    }

    AtomIndicesInside_.reserve(size);

    for (auto it = AtomIndicesBuffer_.beginworker();it != AtomIndicesBuffer_.endworker();it++)
    {
        AtomIndicesInside_.insert(AtomIndicesInside_.end(), it->begin(), it -> end());
    }
}

void DensityField::CalculateInstantaneousInterface()
{
    for (auto it = FieldBuffer_.beginworker();it != FieldBuffer_.endworker();it++)
    {
        it -> zero();
    } 
    // set up the omp buffers
    FieldBuffer_.set_master_object(field_);

    findAtomsIndicesInBoundingBox(); 

    #ifdef MY_DEBUG
    std::cout << "Number of Atoms inside the observation vol is " << AtomIndicesInside_.size() << std::endl;
    #endif 

    auto atomgroup = getAtomGroup(atomGroupName_);
    auto& atoms = atomgroup.getAtoms();

    #pragma omp parallel
    {
        auto& fieldbuf = FieldBuffer_.access_buffer_by_id();

        #pragma omp for 
        for (int i=0;i<AtomIndicesInside_.size();i++)
        {
            int indices = AtomIndicesInside_[i];;
            Real3 correctedPos = bound_box_->PutInBoundingBox(atoms[indices].position);
            index3 Index       = fieldbuf.getClosestGridIndex(correctedPos);


            for(int j=0;j<offsetIndex_.size();j++)
            {
                index3 RealIndex;
                for (int k=0;k<3;k++)
                {
                    RealIndex[k] = Index[k] + offsetIndex_[j][k];
                }

                Real3 latticepos = fieldbuf.getPositionOnGrid(RealIndex[0], RealIndex[1], RealIndex[2]);

                Real3 distance;
                bound_box_->calculateDistance(latticepos, correctedPos, distance);


                Real val = GaussianCoarseGrain::GaussianCoarseGrainFunction(distance, sigma_);
                fieldbuf(RealIndex[0], RealIndex[1], RealIndex[2]) += val;
            }
        }    
    }

    #pragma omp parallel for
    for (int i=0;i<field_.totalSize();i++)
    {
        for (auto it = FieldBuffer_.beginworker();it != FieldBuffer_.endworker();it++)
        {
            field_.accessField()[i] += it->accessField()[i];
        }
    }
}