#include <tools/CommonTypes.h>

#include <cmath>

namespace LinAlg3x3
{
    using Real = CommonTypes::Real;
    using Matrix = CommonTypes::Matrix;
    using Real3= CommonTypes::Real3;

    Real3 CrossProduct(const Real3& v1, const Real3& v2);
    Real DotProduct(const Real3& v1, const Real3& v2);
    Real norm(const Real3& v1);
    Matrix dyad(const Real3& v1, const Real3& v2);
}