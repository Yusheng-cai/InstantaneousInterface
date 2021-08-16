#include "tools/CommonTypes.h"
#include "Curvature.h"
#include "LinAlgTools.h"
#include "Eigen/Core"
#include "Eigen/Dense"

#include <vector>
#include <array>
#include <algorithm>

class CurvatureTensor: public Curvature
{
    public:
        using Matrix = CommonTypes::Matrix;
        using Matrix2= CommonTypes::Matrix2;
        using Real2  = CommonTypes::Real2;
        using Real3  = CommonTypes::Real3;

        CurvatureTensor(CurvatureInput& input); 

        virtual void calculate();
        virtual void printOutput();

        // This is the rotation matrix that rotates vector 1 onto vector 2
        Matrix getRotationMatrix(const Real3& vec1, const Real3& vec2);

        // function that projects curvature from one frame of reference to another
        // project from oldu & oldv to newu & newv
        Real3 projectCurvature(const Real3& oldu, const Real3& oldv, const Real3& refu, const Real3& refv,const Real3& curvature);

    private:
        std::vector<Real3> curvatureTensorPerTriangle_;
        std::vector<Real3> curvatureTensorPerVertex_;

        std::vector<Real2> curvatureVec_;
        std::vector<Real> TotalAreaPerVertex_;
};