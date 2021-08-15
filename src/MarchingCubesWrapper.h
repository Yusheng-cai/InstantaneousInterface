#pragma once
#include "MarchingCubes/TriangulateGlobal.hpp"
#include "Field.h"
#include "tools/CommonTypes.h"
#include "Mesh.h"

#include <vector>
#include <string>
#include <cmath>


class MarchingCubesWrapper
{
    public:
        using Real = CommonTypes::Real;
        using Real3= CommonTypes::Real3;
        using index3= CommonTypes::index3;

        MarchingCubesWrapper(){};
        ~MarchingCubesWrapper(){};

        void calculate(Field& field, Mesh& mesh_, Real isosurfaceVal);
    private:
        TriangulateGlobal marchingCubesalgo_;
};