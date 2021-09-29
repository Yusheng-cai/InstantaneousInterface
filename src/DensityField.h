#pragma once
#include "AtomGroup.h"
#include "tools/InputParser.h"
#include "tools/CommonTypes.h"
#include "SimulationState.h"
#include "BoundingBox.h"
#include "Field.h"
#include "parallel/OpenMP_buffer.h"
#include "GaussianCoarseGrainFunction.h"
#include "tools/GenericFactory.h"
#include "MarchingCubesWrapper.h"
#include "tools/OutputFunction.h"
#include "Curvature.h"


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
};

class DensityField
{
    public:
        using CurvaturePtr = std::unique_ptr<Curvature>;
        using Real = CommonTypes::Real;
        using Real3= CommonTypes::Real3;
        using index3=std::array<int,3>;
        using Range= CommonTypes::Real2;
        using outputfunc = std::function<void(std::string)>;

        DensityField(const DensityFieldInput& input);
        virtual ~DensityField(){};

        virtual void update() {};
        virtual void calculate() = 0;
        virtual void finishCalculate() {};
        virtual void printOutputIfOnStep() {};
        virtual void printFinalOutput();

        void initializeCurvature();

        void findAtomsIndicesInBoundingBox();

        bool isOpen();

        void CalcOffsetIndex();
        void CalculateInstantaneousInterface();

        void addAtomGroup(std::string& name);
        void registerAtomGroupID(std::string& name, int index);

        int getAtomGroupID(std::string& name);
        const AtomGroup& getAtomGroup(std::string& name);
        AtomGroup& accessAtomGroup(std::string& name);

    protected:
        // the atomgroups
        std::vector<const AtomGroup*> AtomGroups_;

        std::map<std::string, int> MapAtomGroupName2Id_;

        std::array<int,3> dimensions_;

        std::string atomGroupName_;

        std::vector<OP::Atom> Allatoms_;

        // The simulation state that keeps track of atom groups and bounding box
        SimulationState& simstate_;

        // sigma used for gaussian smoothing 
        Real sigma_;

        // we cut off n sigmas away
        Real n_ = 2.5;
        Real cutoff_;

        // Name of the bounding box
        std::string boundingboxName_;
        const BoundingBox* bound_box_ = nullptr;

        // Density field related things
        Field field_;

        Range x_range_;
        Range y_range_;
        Range z_range_;

        // offset indices 
        std::vector<index3> offsetIndex_;
        OpenMP::OpenMP_buffer<Field> FieldBuffer_;

        // The isosurface value, usually in units of atom/nm3
        Real isoSurfaceVal_;

        MarchingCubesWrapper MarchingCubes_;
        Mesh mesh_;

        std::vector<int> AtomIndicesInside_;
        OpenMP::OpenMP_buffer<std::vector<int>> AtomIndicesBuffer_;

        // this one manages all the output functions
        Output outputs_;

        // output names 
        std::vector<std::string> OutputNames_;
        std::vector<std::string> OutputFileNames_;

        // ParamterPack 
        ParameterPack& pack_;

        // curvature calculators 
        std::vector<CurvaturePtr> curvatures_;
};

namespace DensityFieldRegistry
{
    using Base = DensityField;
    using Key  = std::string;

    using Factory = GenericFactory<Base, Key, const DensityFieldInput&>;

    template<typename D>
    using registry_ = RegisterInFactory<Base,D, Key, const DensityFieldInput&>;
};