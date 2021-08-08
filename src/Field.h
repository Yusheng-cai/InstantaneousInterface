#pragma once
#include "tools/CommonTypes.h"
#include "tools/Assert.h"

#include <vector>

class Field
{
    public:
        using Real = CommonTypes::Real;
        using Real3= CommonTypes::Real3;
        using Range= CommonTypes::Real2;
        using index3 = CommonTypes::index3;

        Field(){};
        void resize(std::size_t Nx, std::size_t Ny, std::size_t Nz, Range& xrange, Range& yrange, Range& zrange);

        Real& operator()(int i, int j, int k);

        void setXrange(Range& xrange) 
        {   
            x_range_ = xrange;
            Lx_ = x_range_[1] - x_range_[0];
            center_[0] = (x_range_[0] + x_range_[1])*0.5;
        }

        void setYrange(Range& yrange) 
        { 
            y_range_ = yrange;
            Ly_ = y_range_[1] - y_range_[0];
            center_[1] = (y_range_[1] + y_range_[0])*0.5;
        }

        void setZrange(Range& zrange) 
        { 
            z_range_ = zrange; 
            Lz_ = z_range_[1] - z_range_[0];
            center_[2] = (z_range_[1] + z_range_[0])*0.5;
        }

        void fixIndex(index3& index);

        // zero the field
        void zero() {std::fill(field_.begin(), field_.end(),0.0);}

        Real getdx() const {return dx_;}
        Real getdy() const {return dy_;}
        Real getdz() const {return dz_;}


        Real3 getPositionOnGrid(int i, int j, int k);

        // This function assumes that the position has already been PBC corrected to the center of the bb
        index3 getClosestGridIndex(const Real3& position);

    private:
        std::vector<Real> field_;

        int Nx_;
        int Ny_;
        int Nz_;

        Real Lx_ = 0;
        Real Ly_ = 0;
        Real Lz_ = 0;

        Real dx_ = 0.0;
        Real dy_ = 0.0;
        Real dz_ = 0.0;

        Range x_range_ = {{ 0,0 }};
        Range y_range_ = {{ 0,0 }};
        Range z_range_ = {{ 0,0 }};

        Real3 center_;
};