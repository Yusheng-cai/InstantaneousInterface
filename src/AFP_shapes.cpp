#include "AFP_shapes.h"
#include <cmath>

SuperEgg::SuperEgg(const ParameterPack& pack) : AFP_shape(pack){
    pack.ReadNumber("a", ParameterPack::KeyType::Required, a_);
    pack.ReadNumber("b", ParameterPack::KeyType::Required, b_);
    pack.ReadNumber("n", ParameterPack::KeyType::Optional, n_);
    pack.ReadNumber("zmax", ParameterPack::KeyType::Optional, zmax_);
    pack.ReadNumber("a_taper", ParameterPack::KeyType::Optional, a_taper_);
    pack.ReadNumber("b_taper", ParameterPack::KeyType::Optional, b_taper_);
    pack.ReadNumber("a_alpha", ParameterPack::KeyType::Optional, a_alpha_);
    pack.ReadNumber("b_alpha", ParameterPack::KeyType::Optional, b_alpha_);
    pack.ReadArrayNumber("center", ParameterPack::KeyType::Required, center_);
    pack.ReadNumber("offset_height", ParameterPack::KeyType::Optional, offset_height);
}

double SuperEgg::aBulging(double z){
    return a_ * (1 - std::pow((a_taper_ * z / zmax_),2)) + a_alpha_ * std::sin(z / zmax_ * Constants::PI);
}

double SuperEgg::da_dz(double z){
    return a_ * a_taper_ * z / zmax_ * a_taper_ / zmax_ + a_alpha_ * std::cos(z / zmax_ * Constants::PI) * Constants::PI / zmax_;
}


double SuperEgg::bBulging(double z){
    return b_ * (1 - std::pow((b_taper_ * z / zmax_),2)) + b_alpha_ * std::sin(z / zmax_ * Constants::PI);
}

double SuperEgg::db_dz(double z){
    return b_ * (b_taper_ * z / zmax_);
}


SuperEgg::double3 SuperEgg::calculatePos(double u, double v){
    double3 ret;
    ret[2] = zmax_ * std::sin(v);
    ret[0] = aBulging(ret[2]) * Algorithm::sgn(std::cos(u)) * std::pow(std::abs(std::cos(u)), 2.0/(float)n_) + center_[0];
    ret[1] = bBulging(ret[2]) * Algorithm::sgn(std::sin(u)) * std::pow(std::abs(std::sin(u)), 2.0/(float)n_) + center_[1];

    return ret;
}

SuperEgg::double3 SuperEgg::dr_du(double u, double v){
    double3 drdu;
    double dydu, dxdu;
    double z = zmax_ * std::sin(v);

    if (u >= 0){
        dydu = bBulging(z) * 2.0 / (float)n_ * std::pow(std::sin(u), 2.0/(float)n_ - 1) * std::cos(u);
        dxdu = aBulging(z) * 2.0 / (float)n_ * std::pow(std::cos(u), 2.0/(float)n_ - 1) * (-std::sin(u));
    }
    else{
        dydu = bBulging(z) * 2.0 / (float)n_ * std::pow(-std::sin(u), 2.0/(float)n_-1) * (-std::cos(u));
        dxdu = aBulging(z) * 2.0 / (float)n_ * std::pow(-std::cos(u), 2.0/(float)n_-1) * (std::sin(u));
    }

    drdu[0] = dxdu;
    drdu[1] = dydu; 
    drdu[2] = 0.0;

    return drdu;
}

bool SuperEgg::CalculateNormalAndTangent(Real u, Real v, Real3& tangent, Real3& normal, int xdir, int ydir, int zdir){
    // // // initialize the position matrix 
    std::vector<std::vector<double3>> pos(3, std::vector<double3>(3,{0,0,0}));
    // first find the zenithal and azimuthal 
    std::vector<double> ulist = {u - 0.01f, u , u + 0.01f};
    std::vector<double> vlist = {v - 0.01f, v , v + 0.01f};

    // matrix goes like
    // (v-1, u-1), (v-1, u), (v-1, u+1)
    // (v , u-1), (v, u), (v, u+1)
    // (v+1, u-1), (v+1, u), (v+1, u+1)
    
    for (int i=0;i<ulist.size();i++){
        for (int j=0;j<vlist.size();j++){
            pos[j][i] = calculatePos(ulist[i], vlist[j]);
        }
    }


    // these give the dr/du, dr/dv
    Real3 drdu, drdv;
    for (int i=0;i<3;i++){
        drdu[i] = 0.5 * (pos[1][0][i] - pos[1][2][i]);
        drdv[i] = 0.5 * (pos[0][1][i] - pos[2][1][i]);
    }


    Real3 norm = LinAlg3x3::CrossProduct(drdu, drdv);
    LinAlg3x3::normalize(norm);
    normal[xdir] = norm[0];
    normal[ydir] = norm[1];
    normal[zdir] = norm[2];

    Real3 tang = drdv;
    LinAlg3x3::normalize(tang);
    tangent[xdir] = -tang[0];
    tangent[ydir] = -tang[1];
    tangent[zdir] = -tang[2];

    return true;
}




bool SuperEgg::CalculateNormalAndTangent(const double3& point, Real3& tangent, Real3& normal, int xdir, int ydir, int zdir){
    // capturing guess by value, func by reference
    // auto f = [this, point, zdir](Real& x, Real& fx){
    //     fx = this->bBulging(point[zdir] - 1 + 0.864) * Algorithm::sgn(std::sin(x)) * std::pow(std::abs(std::sin(x)), 2.0 / this->getn());

    //     return;
    // };

    // auto f_vec = [this, point](std::vector<Real>& x, std::vector<Real>& fx){
    //     Real z = zmax_ * std::sin(x[0]);
    //     fx[0] = this->aBulging(z) * Algorithm::sgn(std::cos(x[1])) * std::pow(std::abs(std::cos(x[1])), 2.0 / this->getn());
    //     fx[1] = this->bBulging(z) * Algorithm::sgn(std::sin(x[1])) * std::pow(std::abs(std::sin(x[1])), 2.0 / this->getn());

    //     return;
    // };

    // std::function<Real(Real)> f1 = [this, point, zdir, xdir, ydir](Real x){
    //     return this->bBulging(point[zdir]-1+0.864) * Algorithm::sgn(std::sin(x)) * std::pow(std::abs(std::sin(x)), 2.0/this->getn()) - (point[ydir] - this->getcenter()[0]);
    // };


    // first get a good guess of what u and v are
    Real u = std::atan2(point[ydir] - center_[1], point[xdir] - center_[0]);
    if (point[zdir] + offset_height > zmax_){
        return false;
    }
    Real v = std::asin((point[zdir] + offset_height) / zmax_);

    // have a list of u solutions
    std::vector<double> u_solutions;
    std::vector<double> err_list;

    std::vector<double3> pos_s;

    // solve y first 
    FunctorY funcy(point, this, xdir, ydir, zdir);
    Eigen::HybridNonLinearSolver<FunctorY> NL_solver(funcy);
    Eigen::VectorXd y1,y2;
    double3 y1_pos, y2_pos;
    double err_y1, err_y2;
    y1.resize(1), y2.resize(1);
    y1[0] = Algorithm::sgn(u) * Constants::PI  - Algorithm::sgn(u) * 0.1;
    y2[0] = 0 + Algorithm::sgn(u) * 0.1;
    auto infoy = NL_solver.solveNumericalDiff(y1);
    auto infoy2= NL_solver.solveNumericalDiff(y2);

    // keep track of solutions
    u_solutions.push_back(y1[0]);
    u_solutions.push_back(y2[0]);
    y1_pos = calculatePos(y1[0], v), y2_pos = calculatePos(y2[0], v);
    pos_s.push_back(y1_pos); pos_s.push_back(y2_pos);
    err_y1 = std::pow(y1_pos[0] - point[xdir],2) + std::pow(y1_pos[1] - point[ydir],2);
    err_y2 = std::pow(y2_pos[0] - point[xdir],2) + std::pow(y2_pos[1] - point[ydir],2);
    err_list.push_back(err_y1);
    err_list.push_back(err_y2);


    // std::cout << "-----solvey------" << std::endl;
    // std::cout << "status = " << infoy << std::endl;
    // std::cout << "status2 = " << infoy2 << std::endl;
    // std::cout << "Solved u = " << y1[0] << " initial guess = " << u << std::endl;
    // std::cout << "Solved u2= " << y2[0] << std::endl;
    // std::cout << "Solved pos1 = " << y1_pos << std::endl;
    // std::cout << "Solved pos2 = " << y2_pos << std::endl;
    // std::cout << "err_y1 = " << err_y1 << std::endl;
    // std::cout << "err_y2 = " << err_y2 << std::endl;

    // solve x 
    FunctorX funcx(point,this, xdir, ydir, zdir);
    Eigen::HybridNonLinearSolver<FunctorX> NL_solverX(funcx);
    Eigen::VectorXd x1,x2;
    double3 x1_pos, x2_pos;
    double err_x1, err_x2;
    x1.resize(1), x2.resize(1);
    x1[0] = Algorithm::sgn(u) * Constants::PI  - Algorithm::sgn(u) * 0.1;
    x2[0] = 0 + Algorithm::sgn(u) * 0.1;
    auto infox = NL_solverX.solveNumericalDiff(x1);
    auto infox2= NL_solverX.solveNumericalDiff(x2);

    // keep track of solutions
    u_solutions.push_back(x1[0]);
    u_solutions.push_back(x2[0]);
    x1_pos = calculatePos(x1[0],v), x2_pos = calculatePos(x2[0],v);
    pos_s.push_back(x1_pos); pos_s.push_back(x2_pos);
    err_x1 = std::pow(x1_pos[0] - point[xdir],2) + std::pow(x1_pos[1] - point[ydir],2);
    err_x2 = std::pow(x2_pos[0] - point[xdir],2) + std::pow(x2_pos[1] - point[ydir],2);
    err_list.push_back(err_x1);
    err_list.push_back(err_x2);

    // std::cout << "-----solvex------" << std::endl;
    // std::cout << "status = " << infox << std::endl;
    // std::cout << "status2 = " << infox2 << std::endl;
    // std::cout << "Solved u = " << fmod(x1[0], Constants::PI) << " initial guess = " << u << std::endl;
    // std::cout << "Solved u2= " << fmod(x2[0], Constants::PI) << std::endl;
    // std::cout << "Solved pos1 = " << x1_pos << std::endl;
    // std::cout << "Solved pos2 = " << x2_pos << std::endl;
    // std::cout << "err_x1 = " << err_x1 << std::endl;
    // std::cout << "err_x2 = " << err_x2 << std::endl;
    // std::cout << "original pos = " << point << std::endl;

    // set u 
    int index = Algorithm::argmin(err_list);
    u = u_solutions[index];


    // // // initialize the position matrix 
    std::vector<std::vector<double3>> pos(3, std::vector<double3>(3,{0,0,0}));
    // first find the zenithal and azimuthal 
    std::vector<double> ulist = {u - 0.01f, u , u + 0.01f};
    std::vector<double> vlist = {v - 0.01f, v , v + 0.01f};

    // matrix goes like
    // (v-1, u-1), (v-1, u), (v-1, u+1)
    // (v , u-1), (v, u), (v, u+1)
    // (v+1, u-1), (v+1, u), (v+1, u+1)
    
    for (int i=0;i<ulist.size();i++){
        for (int j=0;j<vlist.size();j++){
            pos[j][i] = calculatePos(ulist[i], vlist[j]);
        }
    }


    // these give the dr/du, dr/dv
    Real3 drdu, drdv;
    for (int i=0;i<3;i++){
        drdu[i] = 0.5 * (pos[1][0][i] - pos[1][2][i]);
        drdv[i] = 0.5 * (pos[0][1][i] - pos[2][1][i]);
    }


    Real3 norm = LinAlg3x3::CrossProduct(drdu, drdv);
    LinAlg3x3::normalize(norm);
    normal[xdir] = norm[0];
    normal[ydir] = norm[1];
    normal[zdir] = norm[2];

    Real3 tang = drdv;
    LinAlg3x3::normalize(tang);
    tangent[xdir] = -tang[0];
    tangent[ydir] = -tang[1];
    tangent[zdir] = -tang[2];

    return true;
}