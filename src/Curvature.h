#pragma once
#include "Mesh.h"
#include "tools/CommonTypes.h"
#include "tools/GenericFactory.h"
#include "tools/Assert.h"
#include "tools/InputParser.h"

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

        Curvature(CurvatureInput& input); 

        void close();

        const std::vector<Real>& getCurvature() const {return curvature_;}
        std::vector<Real>& accessCurvature() {return curvature_;}

        // calculates the curvature using a particular method
        virtual void calculate() = 0;
        virtual void printOutput() {};

    protected:
        Mesh& mesh_;

        std::vector<Real> curvature_;

        std::string OutputName_;

        std::ofstream ofs_;
};


namespace CurvatureRegistry
{
    using Base = Curvature;
    using Key  = std::string;

    using Factory = GenericFactory<Base, Key, CurvatureInput&>;

    template<typename D>
    using registry= RegisterInFactory<Base, D, Key, CurvatureInput&>;
};