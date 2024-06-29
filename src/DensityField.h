#pragma once

#ifdef CUDA_ENABLED
#include "GPUArray.cuh"
#include "DensityFieldKernel.h"
#endif

#include "AtomGroup.h"
#include "tools/InputParser.h"
#include "tools/CommonTypes.h"
#include "SimulationState.h"
#include "BoundingBox.h"
#include "Field.h"
#include "parallel/OpenMP_buffer.h"
#include "tools/GenericFactory.h"
#include "marching_cubes.hpp"
#include "tools/OutputFunction.h"
#include "Curvature.h"
#include "Registry.h"
#include "tools/Constants.h"
#include "tools/CommonOperations.h"
#include "surface_nets.h"

#include <vector>
#include <map>
#include <array>
#include <string>
#include <memory>
#include <iostream>
#include <sstream>
#include <functional>

struct DensityFieldInput
{
    SimulationState& simstate_;
    ParameterPack& pack_;
    Registry& reg_;
};

class DensityField
{
    public:
        using CurvaturePtr = std::unique_ptr<Curvature>;
        using Real = CommonTypes::Real;
        using Real3= CommonTypes::Real3;
        using INT3 =std::array<int,3>;
        using INT2 = CommonTypes::index2;
        using Range= CommonTypes::Real2;
        using Meshptr = std::unique_ptr<Mesh>;
        using outputfunc = std::function<void(std::string)>;

        DensityField(const DensityFieldInput& input);
        virtual ~DensityField(){};

        virtual void update() {};
        virtual void calculate() = 0;
        virtual void finishCalculate() {};
        virtual void printOutputIfOnStep() {};
        virtual void printFinalOutput(bool bootstrap=false, int numTimes=0);

        // initialize curvature calculations 
        void initializeCurvature();

        // initialize mesh refinement protocals
        void initializeRefinement();

        // initialize the mesh object
        void initializeMesh();

        // helper function that finds atom indices that lies within the bounding box 
        void findAtomsIndicesInBoundingBox();

        // helper function that calculates the offset index 
        void CalcOffsetIndex();

        // calculate the instantaneous concentration field 
        void CalculateInstantaneousField();

        // calculate the instantaneous concentration field but now is with some property 
        // e.g. q_i * sigma 
        void CalculateInstantaneousFieldProperty(const std::vector<Real>& property);

        // function that adds the atom group to the current object 
        void addAtomGroup(std::string& name);
        void registerAtomGroupID(std::string& name, int index);
        int getAtomGroupID(std::string& name);
        const AtomGroup& getAtomGroup(std::string& name);
        AtomGroup& accessAtomGroup(std::string& name);

        // mesh printing functions
        void printSTL(std::string name);
        void printPLY(std::string name);
        void printBoundaryVertices(std::string name);
        void printnonPBCMesh(std::string name);

        // access the mesh object 
        Mesh& accessMesh() {return *mesh_;}

        // density reset 
        void reset();

        // inline function that coarse grains the concentration
        inline Real GaussianCoarseGrainFunction(const Real3& dx);

    protected:
        // the atomgroups
        std::vector<const AtomGroup*> AtomGroups_;

        std::map<std::string, int> MapAtomGroupName2Id_;

        std::array<int,3> dimensions_;

        std::vector<std::string> atomGroupNames_;

        std::vector<OP::Atom> Allatoms_;

        // The simulation state that keeps track of atom groups and bounding box
        SimulationState& simstate_;

        // registry that keeps track of all the miscellaneous things like curvature
        Registry& reg_;

        // sigma used for gaussian smoothing 
        Real sigma_;
        Real sigmasq_;
        Real prefactor_;
        Real inv_sigmasq2_;

        // we cut off n sigmas away
        Real n_ = 2.5;
        Real cutoff_;

        // whether we do pbc or not for marching cubes 
        bool MCpbc_=false;
        bool useMC_=true;

        // Name of the bounding box
        std::string boundingboxName_;
        const BoundingBox* bound_box_ = nullptr;

        // Density field related things
        Field field_;
        Field property_field_;

        Range x_range_;
        Range y_range_;
        Range z_range_;

        // offset indices 
        std::vector<INT3> offsetIndex_;
        OpenMP::OpenMP_buffer<Field> FieldBuffer_;
        OpenMP::OpenMP_buffer<Field> PropertyFieldBuffer_;

        // The isosurface value, usually in units of atom/nm3
        Real isoSurfaceVal_;

        MarchingCubes MarchingCubes_;
        Surface_nets SurfaceNets_;

        // Mesh stuff 
        Meshptr mesh_;
        std::vector<std::string> MeshOutputs_;
        std::vector<std::string> MeshOutputNames_;
        bool pbcmesh_=false;

        std::vector<std::vector<int>> AtomIndicesInside_;

        // this one manages all the output functions
        Output outputs_;

        // output names 
        std::vector<std::string> OutputNames_;
        std::vector<std::string> OutputFileNames_;

        // ParamterPack 
        ParameterPack& pack_;

        // curvature calculators 
        std::vector<Curvature*> curvatures_;

        // mesh refinement 
        std::vector<MeshRefineStrategy*> refinementstrat_;

        // check whether or not we are cutting a mesh 
        Real3 cut_vec_;
        bool cut_mesh_=false;
        bool cut_below_=false;

        // CUDA stuff --> only defined if cuda is enabled 
        #ifdef CUDA_ENABLED
        GPUArray1d<Real> _atom_property;
        GPUArray2d<Real> _vector_field_neighbors;
        GPUArray2d<int> _vector_field_neighbor_index;
        GPUArray2d<Real> _atom_positions;
        GPUArray3d<Real> _instantaneous_field;
        GPUArray3d<Real> _instantaneous_property_field;
        GPUArray3d<Real> _field;
        GPUArray3d<Real> _property_field;
        GPUArray1d<Real> _box;
        GPUArray1d<Real> _dx;
        GPUArray1d<int> _N;
        GPUArray2d<int> _NeighborIndex;
        int _num_atoms;
        #endif
};

namespace DensityFieldRegistry
{
    using Base = DensityField;
    using Key  = std::string;

    using Factory = GenericFactory<Base, Key, const DensityFieldInput&>;

    template<typename D>
    using registry_ = RegisterInFactory<Base,D, Key, const DensityFieldInput&>;
};