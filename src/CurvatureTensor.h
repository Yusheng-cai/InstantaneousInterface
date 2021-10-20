#include "tools/CommonTypes.h"
#include "Curvature.h"
#include "LinAlgTools.h"
#include "Eigen/Core"
#include "Eigen/Dense"
#include "Graph.h"

#include <vector>
#include <array>
#include <algorithm>
#include <iostream>

class CurvatureTensor: public Curvature
{
    public:
        using Matrix = CommonTypes::Matrix;
        using Matrix2= CommonTypes::Matrix2;
        using Real2  = CommonTypes::Real2;
        using Real3  = CommonTypes::Real3;

        CurvatureTensor(CurvatureInput& input); 
        virtual ~CurvatureTensor(){};

        virtual void calculate();

        virtual void printCurvature(std::string name) override;
        void printPrincipalDirection(std::string name);
        void printFF2(std::string name);
        void printCurvatureDir(std::string name);

        // function that projects curvature from one frame of reference to another
        // project from oldu & oldv to newu & newv
        Real3 projectCurvature(const Real3& oldu, const Real3& oldv, const Real3& refu, const Real3& refv,const Real3& curvature);

        void calculatePrincipalCurvatures();

        // Calculate the curvature in a certain direction as specified by user
        void calculateCurvatureInDir();

    private:
        std::vector<Real3> curvatureTensorPerTriangle_;
        std::vector<Real3> curvatureTensorPerVertex_;
        std::vector<std::array<Real3,2>> PrincipalDirections_;

        std::vector<Real2> curvatureVec_;
        std::vector<Real> TotalAreaPerVertex_;

        std::ofstream curvatureDirOutputofs_;

        std::string curvatureDirOutputName_;

        std::vector<std::vector<Real3>> CurvaturePerVertex_tot;
        std::vector<Real3> curvatureDir_;
        std::vector<std::vector<Real>> curvaturesInDir_;
};