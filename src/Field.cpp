#include "Field.h"

void Field::resize(std::size_t Nx, std::size_t Ny, std::size_t Nz, Range& xrange, Range& yrange, Range& zrange)
{
    ASSERT((Nx > 0), "Nx is not larger than 0 and it is " << Nx << ", the input dimension must be larger than 0.");
    ASSERT((Ny > 0), "Ny is not larger than 0 and it is " << Ny << ", the input dimension must be larger than 0.");
    ASSERT((Nz > 0), "Nz is not larger than 0 and it is " << Nz << ", the input dimension must be larger than 0.");

    Nx_ = Nx;
    Ny_ = Ny;
    Nz_ = Nz;

    setXrange(xrange);
    setYrange(yrange);
    setZrange(zrange);

    field_.resize(Nx*Ny*Nz);

    dx_ = Lx_/Nx_;
    dy_ = Ly_/Ny_;
    dz_ = Lz_/Nz_;
}

Field::Real& Field::operator()(int i, int j, int k)
{
    index3 newIndex = {{i,j,k}}; 
    fixIndex(newIndex);
 
    int index = newIndex[2]*Nx_*Ny_ + newIndex[1]*Nx_ + newIndex[0];

    return field_[index];
}

void Field::fixIndex(index3& index)
{
    ASSERT((index[0] >= -Nx_ && index[0] <= 2*Nx_-1), "x index out of range.");
    ASSERT((index[1] >= -Ny_ && index[1] <= 2*Ny_-1), "y index out of range.");
    ASSERT((index[2] >= -Nz_ && index[2] <= 2*Nz_-1), "z index out of range.");

    if (index[0] < 0 && index[0] >= -Nx_)
    {
        index[0] = index[0] + Nx_;
    }
    else if (index[0] >= Nx_ && index[0] <= 2*Nx_-1)
    {
        index[0] = index[0] - Nx_;
    }
    
    if (index[1] < 0 && index[1] >= -Ny_)
    {
        index[1] = index[1] + Ny_;
    }
    else if (index[1] >= Ny_ && index[1] <= 2*Ny_-1)
    {
        index[1] = index[1] - Ny_;
    }

    if (index[2] < 0 && index[2] >= -Nz_)
    {
        index[2] = index[2] + Nz_;
    }
    else if (index[2] >= Nz_ && index[2] <= 2*Nz_-1)
    {
        index[2] = index[2] - Nz_;
    }
}

Field::Real3 Field::getPositionOnGrid(int i, int j, int k) 
{
    index3 Index = {{i,j,k}};
    fixIndex(Index);

    Real3 ret;

    ret[0] = Index[0]*dx_ + x_range_[0];
    ASSERT((ret[0] <= x_range_[1]), "x is larger than the Bounding Box, x = " << ret[0] << " and bounding box xrange = " << x_range_[1]);

    ret[1] = Index[1]*dy_ + y_range_[0];
    ASSERT((ret[1] <= y_range_[1]), "y is larger than the Bounding Box, y = " << ret[1] << " and bounding box yrange = " << y_range_[1]);

    ret[2] = Index[2]*dz_ + z_range_[0];
    ASSERT((ret[2] <= z_range_[1]), "z is larger than the Bounding Box, z = " << ret[2] << " and bounding box zrange = " << z_range_[1]);

    return ret;
}

Field::index3 Field::getClosestGridIndex(const Real3& position)
{
    index3 ret;
    Real3 PositionNormalized;
    PositionNormalized.fill(0);

    // Normalize the positions by subtracting the least of the xrange, yrange and zrange
    PositionNormalized[0] = position[0] - x_range_[0];
    int xindex = PositionNormalized[0]/dx_;
    ret[0] = xindex;

    PositionNormalized[1] = position[1] - y_range_[0];
    int yindex = PositionNormalized[1]/dy_;
    ret[1] = yindex;

    PositionNormalized[2] = position[2] - z_range_[0];
    int zindex = PositionNormalized[2]/dz_;
    ret[2] = zindex;

    return ret;
}