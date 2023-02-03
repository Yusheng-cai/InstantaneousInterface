#include "tools/CommonTypes.h"
#include "tools/CommonOperations.h"
#include "LinAlgTools.h"
#include <vector>

namespace Contour{
    using Real = CommonTypes::Real;
    using Real3= CommonTypes::Real3;

    std::vector<Real3> ExtendContourLine(const std::vector<Real3>& ContourLine, Real extend_distance);
}