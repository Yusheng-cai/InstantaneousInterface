#pragma once
#include "tools/Assert.h"
#include "tools/InputParser.h"
#include "tools/CommonTypes.h"
#include "tools/GenericFactory.h"
#include "tools/OutputFunction.h"
#include "tools/OutputFunction.h"

#include <vector>
#include <array>


class Mesh;

struct MeshRefineStrategyInput
{
    Mesh& mesh;
    ParameterPack& pack;
};


class MeshRefineStrategy
{
    public:
        using Real = CommonTypes::Real;
        using Real3 = CommonTypes::Real3;

        MeshRefineStrategy(MeshRefineStrategyInput& input);

        virtual ~MeshRefineStrategy(){};

        virtual void refine() = 0;

    protected:
        Mesh& mesh_;
        ParameterPack& pack_;
        Output output_;
};

namespace MeshRefineStrategyFactory
{
    using key = std::string;
    using base = MeshRefineStrategy;

    using factory = GenericFactory<base, key, MeshRefineStrategyInput&>;

    template<typename D>
    using registry_ = RegisterInFactory<base, D, key, MeshRefineStrategyInput&>;
};