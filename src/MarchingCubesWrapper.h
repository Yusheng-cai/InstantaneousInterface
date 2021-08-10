#include "MarchingCubes/TriangulateGlobal.hpp"
#include "Field.h"
#include "tools/CommonTypes.h"

#include <vector>
#include <string>
#include <cmath>

struct vertex
{
    using Real = CommonTypes::Real;
    using Real3= CommonTypes::Real3;

    Real3 position_;
    Real3 normals_;
};

struct triangle
{
    using index3 = CommonTypes::index3;

    index3 triangleindices_;
};

class MarchingCubesWrapper
{
    public:
        using Real = CommonTypes::Real;
        using Real3= CommonTypes::Real3;
        using index3= CommonTypes::index3;

        MarchingCubesWrapper(){};
        ~MarchingCubesWrapper(){};

        void calculate(Field& field, std::vector<vertex>& vertices, std::vector<triangle>& triangles, Real isosurfaceVal);
    private:
        TriangulateGlobal marchingCubesalgo_;
};