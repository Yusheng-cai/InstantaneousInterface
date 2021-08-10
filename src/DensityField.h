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


#include <vector>
#include <map>
#include <array>
#include <string>
#include <memory>
#include <iostream>
#include <sstream>

struct DensityFieldInput
{
    SimulationState& simstate_;
    ParameterPack& pack_;
};

class DensityField
{
    public:
        using Real = CommonTypes::Real;
        using Real3= CommonTypes::Real3;
        using index3=std::array<int,3>;
        using Range= CommonTypes::Real2;

        DensityField(const DensityFieldInput& input);
        ~DensityField(){};

        virtual void update() {};
        virtual void calculate() = 0;
        virtual void finishCalculate() {};
        virtual void printOutputIfOnStep() {};

        bool isOpen();

        void CalcOffsetIndex();

        void addAtomGroup(std::string& name);
        void registerAtomGroupID(std::string& name, int index);

        int getAtomGroupID(std::string& name);
        const AtomGroup& getAtomGroup(std::string& name);
        AtomGroup& accessAtomGroup(std::string& name);

    protected:
        std::vector<const AtomGroup*> AtomGroups_;

        std::map<std::string, int> MapAtomGroupName2Id_;

        std::array<int,3> dimensions_;

        std::string atomGroupName_;

        std::vector<OP::Atom> Allatoms_;

        SimulationState& simstate_;

        Real sigma_;

        std::string output_name_="";
        std::ofstream ofs_;
        std::stringstream ss_;

        // we cut off n sigmas away
        int n_ = 2.5;
        Real cutoff_;

        // Name of the bounding box
        std::string boundingboxName_;
        const BoundingBox* bound_box_ = nullptr;

        // Density field related things
        Field field_;

        Range x_range_;
        Range y_range_;
        Range z_range_;

        std::vector<index3> offsetIndex_;
        OpenMP::OpenMP_buffer<Field> FieldBuffer_;

        Real isoSurfaceVal_;

        MarchingCubesWrapper MarchingCubes_;
        std::vector<triangle> triangles_;
        std::vector<vertex> vertices_;
};

namespace DensityFieldRegistry
{
    using Base = DensityField;
    using Key  = std::string;

    using Factory = GenericFactory<Base, Key, const DensityFieldInput&>;

    template<typename D>
    using registry_ = RegisterInFactory<Base,D, Key, const DensityFieldInput&>;
};