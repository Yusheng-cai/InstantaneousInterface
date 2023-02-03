#include "DensityField.h"

DensityField::DensityField(const DensityFieldInput& input)
:simstate_(input.simstate_), pack_(input.pack_), reg_(input.reg_)
{
    // Read the dimensions of the system, as in (Nx, Ny, Nz)
    input.pack_.ReadArrayNumber("dimensions", ParameterPack::KeyType::Required, dimensions_);

    // Read the atomgroup names of the system
    input.pack_.ReadVectorString("atomgroups", ParameterPack::KeyType::Required, atomGroupNames_);
    for (auto ag : atomGroupNames_){
        addAtomGroup(ag);
    }

    // Read in the sigma value for the gaussian smoothing function
    input.pack_.ReadNumber("sigma", ParameterPack::KeyType::Required, sigma_);
    sigmasq_ = sigma_ * sigma_;
    inv_sigmasq2_ = -1.0 / (2.0 * sigmasq_);
    prefactor_ = std::pow(2.0*Constants::PI*sigmasq_, -1.5);

    // Read in the cut off value (n*sigma_)
    input.pack_.ReadNumber("cutoff", ParameterPack::KeyType::Optional, n_);

    // Read in the isoSurface value
    input.pack_.ReadNumber("isosurfacevalue", ParameterPack::KeyType::Required, isoSurfaceVal_);

    // Read in the pbc 
    input.pack_.Readbool("pbc", ParameterPack::KeyType::Optional, MCpbc_);
    input.pack_.Readbool("useMarchingCubes", ParameterPack::KeyType::Optional, useMC_);

    // Read in the bounding box
    input.pack_.ReadString("boundingbox", ParameterPack::KeyType::Required,boundingboxName_);
    bound_box_ = &simstate_.getBoundingBox(boundingboxName_);
    x_range_ = bound_box_->getXrange();
    y_range_ = bound_box_->getYrange();
    z_range_ = bound_box_->getZrange();

    // read whether or not we are cutting the mesh
    cut_mesh_ = input.pack_.ReadArrayNumber("cut_volume", ParameterPack::KeyType::Optional,cut_vec_);

    // read in the output names 
    input.pack_.ReadVectorString("outputs", ParameterPack::KeyType::Optional, OutputNames_);
    input.pack_.ReadVectorString("outputFiles", ParameterPack::KeyType::Optional, OutputFileNames_);
    ASSERT(( OutputFileNames_.size() == OutputNames_.size()), "The number of outputs does not agree with number of \
    outputs files.");

    // calculate the actual cut off
    cutoff_ = n_*sigma_;

    // resize the field 
    field_.resize(dimensions_[0], dimensions_[1], dimensions_[2], x_range_, y_range_, z_range_);
    property_field_.resize(dimensions_[0],dimensions_[1], dimensions_[2], x_range_, y_range_, z_range_);

    // calculate the offset Index
    CalcOffsetIndex();

    // calculate the field buffer
    for (auto it = FieldBuffer_.beginworker();it != FieldBuffer_.endworker();it++){
        it -> resize(dimensions_[0], dimensions_[1], dimensions_[2], x_range_, y_range_, z_range_);
    }

    for (auto it = PropertyFieldBuffer_.beginworker(); it != PropertyFieldBuffer_.endworker(); it++){
        it -> resize(dimensions_[0], dimensions_[1], dimensions_[2], x_range_, y_range_, z_range_);
    }

    // initialize the mesh object 
    initializeMesh();

    // stores a vector of pointers to the curvature 
    // registry owns the curvature objects 
    initializeCurvature();

    // initialize the refinement process
    initializeRefinement();

    #ifdef CUDA_ENABLED
    // get the maximum number of atoms in each group 
    int size=0;
    int num_neighbors = offsetIndex_.size();
    for (auto ag : atomGroupNames_){
        size += getAtomGroup(ag).getMaxAtoms();
    }

    // resize all the gpu arrays
    _atom_positions.resize(size,3);
    _atom_property.resize(size);
    _vector_field_neighbors.resize(size,num_neighbors);
    _vector_field_neighbor_index.resize(size, num_neighbors);
    _instantaneous_field.resize(dimensions_[0], dimensions_[1], dimensions_[2], 0.0);
    _field.resize(dimensions_[0], dimensions_[1], dimensions_[2], 0.0);
    _property_field.resize(dimensions_[0], dimensions_[1], dimensions_[2], 0.0);
    _instantaneous_property_field.resize(dimensions_[0], dimensions_[1], dimensions_[2], 0.0);
    _NeighborIndex.resize(offsetIndex_.size(),3);
    _box.resize(3);
    _dx.resize(3);
    _N.resize(3);

    // transfer to CUDA
    for (int i=0;i<3;i++){
        _N(i) = dimensions_[i];
    }

    // copy to neighbor index 
    for (int i=0;i<offsetIndex_.size();i++){
        for (int j=0;j<3;j++){
            _NeighborIndex(i,j) = offsetIndex_[i][j];
        }
    }

    _N.memcpyHostToDevice();
    _instantaneous_field.memcpyHostToDevice();
    _field.memcpyHostToDevice();
    _vector_field_neighbors.memcpyHostToDevice();
    _vector_field_neighbor_index.memcpyHostToDevice();
    _NeighborIndex.memcpyHostToDevice();
    #endif

    outputs_.registerOutputFunc("ply", [this](std::string name)-> void {this -> printPLY(name);});
    outputs_.registerOutputFunc("stl", [this](std::string name)-> void {this -> printSTL(name);});
    outputs_.registerOutputFunc("nonpbcMesh", [this](std::string name) -> void {this -> printnonPBCMesh(name);});
}

void DensityField::initializeRefinement()
{
    std::vector<std::string> refinevec;
    pack_.ReadVectorString("refinement", ParameterPack::KeyType::Optional, refinevec);

    for (int i=0;i<refinevec.size();i++){
        auto& r = reg_.getMeshRefineStrat(refinevec[i]);
        refinementstrat_.push_back(&r);
    }
}

void DensityField::initializeCurvature()
{
    std::vector<std::string> cvec;
    pack_.ReadVectorString("curvature", ParameterPack::KeyType::Optional, cvec);

    for (int i=0;i<cvec.size();i++){
        auto& c = reg_.getCurvature(cvec[i]);
        curvatures_.push_back(&c); 
    }
}

void DensityField::reset(){
    // reset the instantaneous field to be zero --> no need to resize 
    #ifdef CUDA_ENABLED
    DensityKernel::FillGPUArray(_field, 0.0);
    #else
    auto& fVec =  field_.accessField();
    std::fill(fVec.begin(), fVec.end(), 0.0);
    #endif
}

void DensityField::printFinalOutput(bool bootstrap, int numTimes){
    if (bootstrap){
        for (int i=0;i<OutputNames_.size();i++){
            std::string name = StringTools::AppendIndexToFileName(OutputFileNames_[i], numTimes);
            outputs_.getOutputFuncByName(OutputNames_[i])(name);
        }
    }
    else{
        for (int i=0;i<OutputNames_.size();i++){
            outputs_.getOutputFuncByName(OutputNames_[i])(OutputFileNames_[i]);
        }
    }

    for (int i=0;i<curvatures_.size();i++){
        curvatures_[i]->printOutput(bootstrap, numTimes);
    }
}

inline DensityField::Real DensityField::GaussianCoarseGrainFunction(const Real3& dx)
{
    Real dotproduct = 0.0;

    for(int i=0;i<3;i++){
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

    for( int i=-Nx_offset; i<=Nx_offset;i++){
        for (int j=-Ny_offset; j<=Ny_offset; j++){
            for (int k=-Nz_offset; k<=Nz_offset;k++){
                INT3 id = {{i,j,k}};
                offsetIndex_.push_back(id);
            }
        }
    }
}

void DensityField::addAtomGroup(std::string& name){
    int index = AtomGroups_.size();

    registerAtomGroupID(name, index);

    AtomGroups_.push_back(&simstate_.getAtomGroup(name));
}

void DensityField::registerAtomGroupID(std::string& name, int index){
    auto it = MapAtomGroupName2Id_.find(name);

    ASSERT(( it == MapAtomGroupName2Id_.end()), "The name for atomgroup " << name << " is already registered.");

    MapAtomGroupName2Id_.insert(std::make_pair(name, index));
}

int DensityField::getAtomGroupID(std::string& name){
    auto it = MapAtomGroupName2Id_.find(name);

    ASSERT(( it != MapAtomGroupName2Id_.end()), "The name for atomgroup " << name << " is already registered.");

    return it -> second;
}

const AtomGroup& DensityField::getAtomGroup(std::string& name){
    int ID = getAtomGroupID(name);

    return *AtomGroups_[ID];
}

AtomGroup& DensityField::accessAtomGroup(std::string& name){
    int ID = getAtomGroupID(name);

    return const_cast<AtomGroup&>(*AtomGroups_[ID]);
}

void DensityField::initializeMesh(){
    // check if we want to output a pbc mesh 
    pack_.Readbool("PBCMesh", ParameterPack::KeyType::Optional, MCpbc_);

    // initialize the mesh object 
    mesh_ = Meshptr(new Mesh());
    if (MCpbc_){
        mesh_->setBoxLength(bound_box_->getSides());
    }
}

void DensityField::findAtomsIndicesInBoundingBox(){
    AtomIndicesInside_.clear();
    for (auto ag : atomGroupNames_){
        const auto& atomgroup = getAtomGroup(ag);
        auto& atoms = atomgroup.getAtoms();

        std::vector<int> indices_ag;

        #pragma omp parallel
        {
            std::vector<int> indices_local; 

            #pragma omp for
            for(int i=0;i<atoms.size();i++){
                if (bound_box_ -> isInside(atoms[i].position)){
                    indices_local.push_back(i);
                }
            }

            #pragma omp critical
            {
                indices_ag.insert(indices_ag.end(), indices_local.begin(), indices_local.end());
            }
        }

        AtomIndicesInside_.push_back(indices_ag);
    }

    #ifdef CUDA_ENABLED
    _num_atoms=0;
    for (int i=0;i<AtomIndicesInside_.size();i++){
        const auto& atomgroup = getAtomGroup(atomGroupNames_[i]);
        auto& atoms = atomgroup.getAtoms();
        for (int j=0;j<AtomIndicesInside_[i].size();j++){
            int ind = AtomIndicesInside_[i][j];
            Real3 correctedPos = bound_box_->PutInBoundingBox(atoms[ind].position);
            for (int k=0;k<3;k++){
                _atom_positions(_num_atoms,k) = correctedPos[k];
            }
            _num_atoms++;
        }
    }

    #endif
}

void DensityField::CalculateInstantaneousFieldProperty(const std::vector<Real>& property){
    findAtomsIndicesInBoundingBox();

    #ifdef CUDA_ENABLED
    ASSERT((property.size() <= _atom_property.getTotalSize()), "atom property vector too small, property size = " \
     << property.size() << " while atom property is resize to " << _atom_property.getTotalSize());
    for (int i=0;i<property.size();i++){
        _atom_property(i) = property[i];
    }

    for (int i=0;i<3;i++){
        _box(i) = simstate_.getSimulationBox().getSides()[i];
        _dx(i)  = _box(i) / _N(i);
    }

    _box.memcpyHostToDevice();
    _dx.memcpyHostToDevice();
    _atom_positions.memcpyHostToDevice();
    _atom_property.memcpyHostToDevice();
    DensityKernel::CalculateInstantaneousFieldProperty(_vector_field_neighbors, _atom_positions, _atom_property, _vector_field_neighbor_index,\
                                                inv_sigmasq2_, prefactor_, _box, \
                                               _dx, _N, _instantaneous_field, _field,_NeighborIndex, _num_atoms);
    #else
    // the master object will not be zero'd
    for (auto it = FieldBuffer_.beginworker();it != FieldBuffer_.endworker();it++){
        it -> zero();
    } 
    for (auto it = PropertyFieldBuffer_.beginworker();it != PropertyFieldBuffer_.endworker();it++){
        it -> zero();
    } 

    // set up the omp buffers
    FieldBuffer_.set_master_object(field_);
    PropertyFieldBuffer_.set_master_object(property_field_);

    // start calculating InstantaneousInterface
    for (int i=0;i<atomGroupNames_.size();i++)
    {
        std::string ag = atomGroupNames_[i];
        const auto& atomgroup = getAtomGroup(ag);
        auto& atoms = atomgroup.getAtoms();

        #pragma omp parallel
        {
            auto& fieldbuf = FieldBuffer_.access_buffer_by_id();
            auto& propertybuf = PropertyFieldBuffer_.access_buffer_by_id();

            #pragma omp for 
            for (int j=0;j<AtomIndicesInside_[i].size();j++){
                int indices        = AtomIndicesInside_[i][j];;
                Real3 correctedPos = bound_box_->PutInBoundingBox(atoms[indices].position);
                Real p       = property[indices];
                INT3  Index        = fieldbuf.getClosestGridIndex(correctedPos);

                // Fix the index --> to correct for periodic boundary conditions  
                fieldbuf.fixIndex(Index);

                // iterate over all the indices and calculate instantaneousinterface
                for(int k=0;k<offsetIndex_.size();k++){
                    // offset index 
                    INT3 RealIndex = Index + offsetIndex_[k];

                    // plate the positions on grid
                    Real3 latticepos = fieldbuf.getPositionOnGrid(RealIndex[0], RealIndex[1], RealIndex[2]);

                    // calculate the distance 
                    Real3 distance;
                    bound_box_->calculateDistance(latticepos, correctedPos, distance);

                    Real gauss = GaussianCoarseGrainFunction(distance);
                    Real val = p * gauss;
                    propertybuf(RealIndex[0], RealIndex[1], RealIndex[2]) += val;
                    fieldbuf(RealIndex[0], RealIndex[1], RealIndex[2]) += gauss;
                }
            }    
        }

        #pragma omp parallel for
        for (int i=0;i<field_.totalSize();i++)
        {
            for (auto it = FieldBuffer_.beginworker();it != FieldBuffer_.endworker();it++){
                field_.accessField()[i] += it->accessField()[i];
            }
        }

        #pragma omp parallel for
        for (int i=0;i<property_field_.totalSize();i++){
            for (auto it = PropertyFieldBuffer_.beginworker(); it != PropertyFieldBuffer_.endworker(); it++){
                property_field_.accessField()[i] += it -> accessField()[i];
            }
        }
    }
    #endif
}

void DensityField::CalculateInstantaneousField(){
    // find all the atom groups indices that are in the bounding box
    findAtomsIndicesInBoundingBox(); 

    #ifdef CUDA_ENABLED
    // get the simulation box 
    for (int i=0;i<3;i++){
        _box(i) = simstate_.getSimulationBox().getSides()[i];
        _dx(i)  = _box(i) / _N(i);
    }

    _box.memcpyHostToDevice();
    _dx.memcpyHostToDevice();
    _atom_positions.memcpyHostToDevice();
    DensityKernel::CalculateInstantaneousField(_vector_field_neighbors, _atom_positions, _vector_field_neighbor_index,\
                                                inv_sigmasq2_, prefactor_, _box, \
                                               _dx, _N, _instantaneous_field, _field,_NeighborIndex, _num_atoms);
    #else
    // the master object will not be zero'd
    for (auto it = FieldBuffer_.beginworker();it != FieldBuffer_.endworker();it++){
        it -> zero();
    } 

    // set up the omp buffers
    FieldBuffer_.set_master_object(field_);

    // start calculating InstantaneousInterface
    for (int i=0;i<atomGroupNames_.size();i++)
    {
        std::string ag = atomGroupNames_[i];
        const auto& atomgroup = getAtomGroup(ag);
        auto& atoms = atomgroup.getAtoms();

        #pragma omp parallel
        {
            auto& fieldbuf = FieldBuffer_.access_buffer_by_id();

            #pragma omp for 
            for (int j=0;j<AtomIndicesInside_[i].size();j++){
                int indices        = AtomIndicesInside_[i][j];;
                Real3 correctedPos = bound_box_->PutInBoundingBox(atoms[indices].position);
                INT3  Index        = fieldbuf.getClosestGridIndex(correctedPos);

                // Fix the index --> to correct for periodic boundary conditions  
                fieldbuf.fixIndex(Index);

                // iterate over all the indices and calculate instantaneousinterface
                for(int k=0;k<offsetIndex_.size();k++){
                    // offset index 
                    INT3 RealIndex = Index + offsetIndex_[k];

                    // plate the positions on grid
                    Real3 latticepos = fieldbuf.getPositionOnGrid(RealIndex[0], RealIndex[1], RealIndex[2]);

                    // calculate the distance 
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
            for (auto it = FieldBuffer_.beginworker();it != FieldBuffer_.endworker();it++){
                field_.accessField()[i] += it->accessField()[i];
            }
        }
    }
    #endif
}

void DensityField::printSTL(std::string name){
    MeshTools::writeSTL(name, *mesh_);
}

void DensityField::printPLY(std::string name){
    MeshTools::writePLY(name, *mesh_);
}

void DensityField::printBoundaryVertices(std::string name){
    std::ofstream ofs;
    ofs.open(name);

    std::vector<bool> boundaryIndicator;
    std::vector<std::vector<INT2>> MapVertexToEdges;
    std::map<INT2, std::vector<int>> MapEdgeToFace;

    MeshTools::MapEdgeToFace(*mesh_, MapEdgeToFace, MapVertexToEdges);
    MeshTools::CalculateBoundaryVertices(*mesh_, MapEdgeToFace, boundaryIndicator);
    
    const auto& v = mesh_->getvertices();

    ofs << "# Index Vx Vy Vz Nx Ny Nz" << "\n";
    for (int i=0;i<v.size();i++){
        if (MeshTools::IsBoundary(i, boundaryIndicator)){
            ofs << i << " " << v[i].position_[0] << " " << v[i].position_[1] << " " << v[i].position_[2] << \
            " " << v[i].normals_[0] << " " << v[i].normals_[1] << " " << v[i].normals_[2] << "\n";
        }
    }

    ofs.close();
}

void DensityField::printnonPBCMesh(std::string name){
    MeshTools::writeNonPBCMesh(name, *mesh_);
}