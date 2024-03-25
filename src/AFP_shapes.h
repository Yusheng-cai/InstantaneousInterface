#pragma once

#include "tools/CommonTypes.h"
#include "tools/InputParser.h"
#include "tools/Constants.h"
#include "LinAlgTools.h"
#include "tools/Algorithm.h"
#include "fsolve_wrapper.hpp"
#include "Eigen/Core"


#include <vector>

class AFP_shape {
    using Real3 = CommonTypes::Real3;
    using Real  = CommonTypes::Real;
    using double3= CommonTypes::double3;

    public:
        AFP_shape(const ParameterPack& pack){};

        virtual void CalculateNormalAndTangent(const double3& point, Real3& tangent, Real3& normal, int xdir=1, int ydir=2, int zdir=0)=0; 
};

class SuperEgg : public AFP_shape{
    public:
        using Real3 = CommonTypes::Real3;
        using Real  = CommonTypes::Real;
        using Real2 = CommonTypes::Real2;
        using double3 = CommonTypes::double3;
        using double2= std::array<double,2>;

        SuperEgg(const ParameterPack& pack);

        double aBulging(double z);
        double bBulging(double z);
        double da_dz(double z);
        double db_dz(double z);

        double3 calculatePos(double u, double v);
        virtual void CalculateNormalAndTangent(const double3& point, Real3& tangent, Real3& normal, int xdir=1, int ydir=2, int zdir=0) override;

        double getn() {return (float)n_;}
        double2 getcenter() {return center_;}

        double3 dr_du(double u, double v);

        struct FunctorY{
            double3 point_;
            SuperEgg* se_;
            FunctorY(double3 point, SuperEgg* se) : point_(point), se_(se) {};

            int operator()(Eigen::VectorXd &b, Eigen::VectorXd &fvec){
                fvec[0] = se_->bBulging(point_[0] + se_->offset_height) * Algorithm::sgn(std::sin(b[0])) * std::pow(std::abs(std::sin(b[0])), 2.0/se_->getn()) - (point_[2] - se_->getcenter()[1]);

                return 0;
            }
            
            int operator()(Eigen::VectorXd &b, Eigen::VectorXd &fvec) const{
                fvec[0] = se_->bBulging(point_[0] + se_->offset_height) * Algorithm::sgn(std::sin(b[0])) * std::pow(std::abs(std::sin(b[0])), 2.0/se_->getn()) - (point_[2] - se_->getcenter()[1]);

                return 0;
            }
        };

        struct FunctorX{
            double3 point_;
            SuperEgg* se_;
            FunctorX(double3 point, SuperEgg* se) : point_(point), se_(se) {};

            int operator()(Eigen::VectorXd &b, Eigen::VectorXd &fvec){
                fvec[0] = se_->aBulging(point_[0] + se_->offset_height) * Algorithm::sgn(std::cos(b[0])) * std::pow(std::abs(std::cos(b[0])), 2.0/se_->getn()) - (point_[1] - se_->getcenter()[0]);

                return 0;
            }
            
            int operator()(Eigen::VectorXd &b, Eigen::VectorXd &fvec) const{
                fvec[0] = se_->aBulging(point_[0] + se_->offset_height) * Algorithm::sgn(std::cos(b[0])) * std::pow(std::abs(std::cos(b[0])), 2.0/se_->getn()) - (point_[1] - se_->getcenter()[0]);

                return 0;
            }
        };

    private:
        double a_;
        double b_;
        int n_=4;
        double zmax_=1.5;
        double a_taper_=1.0;
        double b_taper_=0.5;
        double a_alpha_=0.0;
        double b_alpha_=0.0;
        double offset_height=0.0;
        double2 center_;
};
