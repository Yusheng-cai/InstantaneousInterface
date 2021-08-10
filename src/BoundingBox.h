#pragma once
#include "tools/InputParser.h"
#include "tools/CommonTypes.h"
#include "tools/Assert.h"

#include <vector>
#include <array>
#include <string>

class SimulationState;

struct BoundingBoxInput
{
    ParameterPack& pack_;
    SimulationState& simstate_;
};

class BoundingBox
{
    public:
        using Real = CommonTypes::Real;
        using Real3= CommonTypes::Real3;
        using Range= std::array<Real,2>;

        BoundingBox(const BoundingBoxInput& input);

        const Range& getXrange() const {return x_range_;} 
        const Range& getYrange() const {return y_range_;}
        const Range& getZrange() const {return z_range_;}

        bool isInside(const Real3& position) const;

        void calculateDistance(const Real3& x1, const Real3& x2, Real3& distance) const;

        // A function that shifts the atom with respect to the center of the bounding box
        Real3 PutInBoundingBox(const Real3& position) const;

        std::string getName() const {return name_;}

    
    private:
        Real Lx_;
        Real Ly_;
        Real Lz_;

        Real halfLx_;
        Real halfLy_;
        Real halfLz_;

        Range x_range_;
        Range y_range_;
        Range z_range_;

        SimulationBox& simBox_;
        SimulationBox BoundBox_;

        Real3 center_;

        std::string name_;

        bool calculateWithBoundBox_ = false;
        std::vector<int> ignorePBCdims_;
};