#pragma once

#include "tools/CommonTypes.h"
#include "tools/Assert.h"

#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>

// concentration field 
class Field
{
    public:
        using Real = CommonTypes::Real;
        using Real3= CommonTypes::Real3;
        using Range= CommonTypes::Real2;
        using index3 = CommonTypes::index3;
        using Range3= Real3;

        Field(){};
        Field(std::size_t Nx, std::size_t Ny, std::size_t Nz, Range& xrange, Range& yrange, Range& zrange);
        void resize(std::size_t Nx, std::size_t Ny, std::size_t Nz, Range& xrange, Range& yrange, Range& zrange);
        void clearField();

        void pushback(Real val);

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
        index3 size() const;
        int totalSize() const {return Nx_*Ny_*Nz_;}
        Real* data() {return field_.data();}

        std::vector<Real>& accessField() {return field_;}


        Real getdx() const {return dx_;}
        Real getdy() const {return dy_;}
        Real getdz() const {return dz_;}
        Real3 getLength() const
        {
            Real3 length_;
            length_[0] = Lx_;
            length_[1] = Ly_;
            length_[2] = Lz_;

            return length_;
        }

        Real3 getSpacing() const 
        {
            Real3 space;
            space = {{ dx_, dy_, dz_}};

            return space;
        }

        index3 getN() const 
        {
            index3 index;
            index = {{Nx_, Ny_, Nz_}};

            return index;
        }

        const Range& getXrange() const {return x_range_;}
        const Range& getYrange() const {return y_range_;}
        const Range& getZrange() const {return z_range_;}

        const Real3 getMinRange() const 
        {
            Real3 minRange = {{ x_range_[0], y_range_[0], z_range_[0]}};

            return minRange;
        }

        const Real3 getMaxRange() const 
        {
            Real3 maxRange = {{x_range_[1], y_range_[1], z_range_[1]}};

            return maxRange;
        }


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

namespace FieldTools
{
    using Real = CommonTypes::Real;

    // assume field has already been initialized, as in the xrange, yrange zrange and Nx, Ny Nz has been set 
    void FieldReader(std::string& FileName, Field& f);
};