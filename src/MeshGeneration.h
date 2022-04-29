#pragma once 

#include "tools/Assert.h"
#include "tools/CommonTypes.h"
#include "tools/GenericFactory.h"
#include "tools/InputParser.h"
#include "tools/OutputFunction.h"

#include <iostream>
#include <vector>
#include <array>


class MeshGeneration
{
    public:
        using Real  = CommonTypes::Real;
        using Real3 = CommonTypes::Real3;
        using Real2 = CommonTypes::Real2;
        using INT2  = CommonTypes::index2;
        using INT3  = CommonTypes::index3;

        MeshGeneration(ParameterPack& pack);
        virtual void generate() = 0;
        virtual void printOutputs();

    protected:
        ParameterPack& pack_;
        Output output_;

        std::vector<std::string> OutputNames_;
        std::vector<std::string> OutputFileNames_;

        // neighbor Indices
        std::vector<std::vector<int>> neighborIndices_;
        std::map<INT2, std::vector<int>> MapEdgeToFace_;
        std::vector<std::vector<INT2>> MapVertexToEdge_;
        std::vector<bool> boundaryIndicator_;
};


namespace MeshGenerationFactory
{
    using key = std::string;
    using base = MeshGeneration;
    using factory = GenericFactory<base, key, ParameterPack&>;

    template<typename D>
    using registry_ = RegisterInFactory<base, D, key, ParameterPack&>;
};