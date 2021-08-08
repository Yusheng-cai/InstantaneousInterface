#include "tools/Assert.h"
#include "tools/InputParser.h"
#include "tools/CommonTypes.h"
#include "tools/Constants.h"

#include <cmath>
#include <vector>

namespace GaussianCoarseGrain
{
    using Real = CommonTypes::Real;
    using Real3 = CommonTypes::Real3;

    inline Real GaussianCoarseGrainFunction(const Real3& dx, Real sigma);
}

inline GaussianCoarseGrain::Real GaussianCoarseGrain::GaussianCoarseGrainFunction(const Real3& dx, Real sigma)
{
    Real sigma2 = sigma*sigma;
    Real factor = std::pow(2*Constants::PI*sigma2, -1.5);

    Real dotproduct = 0.0;

    for(int i=0;i<3;i++)
    {
        dotproduct += dx[i] * dx[i];
    }

    return factor*std::exp(-dotproduct/(2*sigma2));
}