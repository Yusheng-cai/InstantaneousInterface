#include "SimulationState.h"

BoundingBox::BoundingBox(const BoundingBoxInput& input)
:simBox_(input.simstate_.getSimulationBox())
{
    input.pack_.ReadString("name", ParameterPack::KeyType::Required, name_);
    input.pack_.ReadArrayNumber("xrange", ParameterPack::KeyType::Required, x_range_);
    input.pack_.ReadArrayNumber("yrange", ParameterPack::KeyType::Required, y_range_);
    input.pack_.ReadArrayNumber("zrange", ParameterPack::KeyType::Required, z_range_);

    ASSERT((x_range_[1] > x_range_[0]), "The range needs to be provided with smaller value at first and they can not be the same value.");
    ASSERT((y_range_[1] > y_range_[0]), "The range needs to be provided with smaller value at first and they can not be the same value.");
    ASSERT((z_range_[1] > z_range_[0]), "The range needs to be provided with smaller value at first and they can not be the same value.");

    Lx_ = x_range_[1] - x_range_[0];
    Ly_ = y_range_[1] - y_range_[0];
    Lz_ = z_range_[1] - z_range_[0];

    halfLx_ = 0.5*Lx_;
    halfLy_ = 0.5*Ly_;
    halfLz_ = 0.5*Lz_;

    center_[0] = 0.5*(x_range_[1] - x_range_[0]) + x_range_[0];
    center_[1] = 0.5*(y_range_[1] - y_range_[0]) + y_range_[0];
    center_[2] = 0.5*(z_range_[1] - z_range_[0]) + z_range_[0];

    sides_ = {{Lx_, Ly_, Lz_}};
}

BoundingBox::Real3 BoundingBox::PutInBoundingBox(const Real3& position) const
{    
    Real3 ret;
    ret.fill(0);

    Real3 shift = simBox_.calculateShift(position, center_);

    for (int i=0;i<3;i++)
    {
        ret[i]  = position[i] + shift[i];
    }

    return ret;
}

bool BoundingBox::isInside(const Real3& position) const
{
    Real3 distance;
    Real sq_dist;

    simBox_.calculateDistance(position, center_, distance, sq_dist);

    if (std::abs(distance[0]) > halfLx_)
    {
        return false;
    }

    if (std::abs(distance[1]) > halfLy_)
    {
        return false;
    }

    if (std::abs(distance[2]) > halfLz_)
    {
        return false;
    }

    return true;
}

void BoundingBox::calculateDistance(const Real3& x1, const Real3& x2, Real3& distance) const
{
    Real sq_dist;
 
    simBox_.calculateDistance(x1, x2, distance, sq_dist);
}