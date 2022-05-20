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
        virtual void calculate(Mesh& mesh) = 0;
        void CalculateFaceCurvature(Mesh& mesh, const std::vector<Real>& VertexCurvature, std::vector<Real>& FaceCurvature);

        virtual void printOutput();
        virtual void printCurvature(std::string name);
        virtual void printPrincipalDir(std::string name);
        virtual void printFaceCurvature(std::string name);

        std::string getName() {return name_;}

        void initialize(Mesh& mesh);

        const std::vector<Real>& getAvgCurvaturePerVertex() {return avgCurvaturePerVertex_;}
        const std::vector<Real>& getGaussCurvaturePerVertex() {return GaussCurvaturePerVertex_;}

    protected:
        ParameterPack& pack_;

        std::string name_="c";

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

        // the curvature for each triangle --> k1, k2 , avgCurvature
        std::vector<Real> FaceCurvature_;
};


namespace CurvatureRegistry
{
    using Base = Curvature;
    using Key  = std::string;

    using Factory = GenericFactory<Base, Key, CurvatureInput&>;

    template<typename D>
    using registry= RegisterInFactory<Base, D, Key, CurvatureInput&>;
};