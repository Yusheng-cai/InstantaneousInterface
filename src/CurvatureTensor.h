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

        CurvatureTensor(CurvatureInput& input); 

        virtual void calculate();
        virtual void printOutput();
    private:
        std::vector<Real3> curvatureTensor_;

        std::vector<Real2> curvatureVec_;
};