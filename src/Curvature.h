#pragma once
#include "Mesh.h"
#include "tools/CommonTypes.h"
#include "tools/GenericFactory.h"
#include "tools/Assert.h"
#include "tools/OutputFunction.h"
#include "tools/InputParser.h"
#include "Eigen/Dense"
#include "Eigen/Core"
#include "happly.h"

#include <vector>
#include <array>
#include <string>
#include <iostream>

struct CurvatureInput
{
    ParameterPack& pack;
    Mesh& mesh_;
};

class Curvature
{
    public:
        using Real = CommonTypes::Real;
        using Real3= CommonTypes::Real3;
        using Real2 = CommonTypes::Real2;

        Curvature(CurvatureInput& input); 
        virtual ~Curvature(){};

        // calculates the curvature using a particular method
        virtual void calculate() = 0;
        virtual void printOutput();
        virtual void printCurvature(std::string name);
        virtual void printPLYlibr(std::string name);
        virtual void printPrincipalDir(std::string name);

    protected:
        Mesh& mesh_;

        // The output files
        Output outputs_;

        std::vector<Real> avgCurvaturePerVertex_;
        std::vector<Real> GaussCurvaturePerVertex_;
        std::vector<Real2> CurvaturePerVertex_;

        std::vector<std::string> OutputNames_;
        std::vector<std::string> OutputFileNames_;

        // The principal directions used for some curvature calculations
        std::vector<Real3> principalDir1_;
        std::vector<Real3> principalDir2_;
};


namespace CurvatureRegistry
{
    using Base = Curvature;
    using Key  = std::string;

    using Factory = GenericFactory<Base, Key, CurvatureInput&>;

    template<typename D>
    using registry= RegisterInFactory<Base, D, Key, CurvatureInput&>;
};