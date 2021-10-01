#pragma once
#include "tools/Assert.h"
#include "tools/InputParser.h"
#include "tools/CommonTypes.h"

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

        virtual void refine() = 0;

    protected:
        Mesh& mesh_;
        ParameterPack& pack_;
};