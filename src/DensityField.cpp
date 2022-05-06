#include "DensityField.h"

DensityField::DensityField(const DensityFieldInput& input)
:simstate_(input.simstate_), pack_(input.pack_), reg_(input.reg_)
{
    // Read the dimensions of the system, as in (Nx, Ny, Nz)
    input.pack_.ReadArrayNumber("dimensions", ParameterPack::KeyType::Required, dimensions_);

    // Read the atomgroup names of the system
    input.pack_.ReadString("atomgroup", ParameterPack::KeyType::Required, atomGroupName_);
    addAtomGroup(atomGroupName_);

    // Read in the sigma value for the gaussian smoothing function
    input.pack_.ReadNumber("sigma", ParameterPack::KeyType::Required, sigma_);
    sigmasq_ = sigma_ * sigma_;
    prefactor_ = std::pow(2*Constants::PI*sigmasq_, -1.5);

    // Read in the cut off value (n*sigma_)
    input.pack_.ReadNumber("cutoff", ParameterPack::KeyType::Optional, n_);

    // Read in the isoSurface value
    input.pack_.ReadNumber("isosurfacevalue", ParameterPack::KeyType::Required, isoSurfaceVal_);

    // Read in the pbc 
    input.pack_.Readbool("pbc", ParameterPack::KeyType::Optional, MCpbc_);

    // Read in the bounding box
    input.pack_.ReadString("boundingbox", ParameterPack::KeyType::Required,boundingboxName_);
    bound_box_ = &simstate_.getBoundingBox(boundingboxName_);
    x_range_ = bound_box_->getXrange();
    y_range_ = bound_box_->getYrange();
    z_range_ = bound_box_->getZrange();

    // read in the output names 
    input.pack_.ReadVectorString("outputs", ParameterPack::KeyType::Optional, OutputNames_);
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

    // stores a vector of pointers to the curvature 
    // registry owns the curvature objects 
    initializeCurvature();

    // initialize the refinement process
    initializeRefinement();
}

void DensityField::initializeRefinement()
{
    std::vector<std::string> refinevec;
    pack_.ReadVectorString("refinement", ParameterPack::KeyType::Optional, refinevec);

    for (int i=0;i<refinevec.size();i++)
    {
        auto& r = reg_.getMeshRefineStrat(refinevec[i]);
        refinementstrat_.push_back(&r);
    }
}

void DensityField::initializeCurvature()
{
    std::vector<std::string> cvec;
    pack_.ReadVectorString("curvature", ParameterPack::KeyType::Optional, cvec);

    for (int i=0;i<cvec.size();i++)
    {
        auto& c = reg_.getCurvature(cvec[i]);
        curvatures_.push_back(&c); 
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

inline DensityField::Real DensityField::GaussianCoarseGrainFunction(const Real3& dx)
{
    Real dotproduct = 0.0;

    for(int i=0;i<3;i++)
    {
        dotproduct += dx[i] * dx[i];
    }

    return prefactor_*std::exp(-dotproduct/(2*sigmasq_));
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
                INT3 id = {{i,j,k}};
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

void DensityField::CalculateInstantaneousField()
{
    // the master object will not be zero'd
    for (auto it = FieldBuffer_.beginworker();it != FieldBuffer_.endworker();it++)
    {
        it -> zero();
    } 

    // set up the omp buffers
    FieldBuffer_.set_master_object(field_);

    findAtomsIndicesInBoundingBox(); 

    auto atomgroup = getAtomGroup(atomGroupName_);
    auto& atoms = atomgroup.getAtoms();

    #pragma omp parallel
    {
        auto& fieldbuf = FieldBuffer_.access_buffer_by_id();

        #pragma omp for 
        for (int i=0;i<AtomIndicesInside_.size();i++)
        {
            int indices        = AtomIndicesInside_[i];;
            Real3 correctedPos = bound_box_->PutInBoundingBox(atoms[indices].position);
            INT3  Index        = fieldbuf.getClosestGridIndex(correctedPos);

            // Fix the index 
            fieldbuf.fixIndex(Index);

            for(int j=0;j<offsetIndex_.size();j++)
            {
                INT3 RealIndex;
                for (int k=0;k<3;k++)
                {
                    RealIndex[k] = Index[k] + offsetIndex_[j][k];
                }

                Real3 latticepos = fieldbuf.getPositionOnGrid(RealIndex[0], RealIndex[1], RealIndex[2]);

                Real3 distance;
                bound_box_->calculateDistance(latticepos, correctedPos, distance);

                Real val = GaussianCoarseGrainFunction(distance);
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