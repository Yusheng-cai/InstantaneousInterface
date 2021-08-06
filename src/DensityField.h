#pragma once
#include "AtomGroup.h"
#include "tools/InputParser.h"
#include "tools/CommonTypes.h"
#include "SimulationState.h"
#include "BoundingBox.h"
#include "Field.h"

#include <vector>
#include <array>
#include <string>
#include <memory>

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

        void addAtomGroup();

    private:
        std::vector<const AtomGroup*> AtomGroups_;

        std::array<int,3> dimensions_;

        std::vector<std::string> atomGroupNames_;

        SimulationState& simstate_;

        Real sigma_;

        // we cut off n sigmas away
        int n_ = 2.5;

        Real cutoff_;

        std::string boundingboxName_;
        const BoundingBox* bound_box_ = nullptr;

        Field field_;

        Range x_range_;
        Range y_range_;
        Range z_range_;
};