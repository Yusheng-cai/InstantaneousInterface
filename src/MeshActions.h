#pragma once

#include <vector>
#include <string>
#include <array>
#include <memory>

#include "tools/Assert.h"
#include "tools/CommandLineArguments.h"
#include "Mesh.h"
#include "tools/CommonTypes.h"
#include "Curvature.h"
#include "tools/InputParser.h"
#include "MeshRefineStrategy.h"
#include "CurvatureCurveFit.h"

namespace MeshActions
{
    using INT3 = std::array<int,3>;
    using Real = CommonTypes::Real;
    using Real3= CommonTypes::Real3;
    using curveptr = std::unique_ptr<Curvature>;
    using refineptr= std::unique_ptr<MeshRefineStrategy>;

    void TranslateMesh(CommandLineArguments& cmd);
    void CurveFit(CommandLineArguments& cmd);
    void JetFit(CommandLineArguments& cmd);
    void CurvatureFlow(CommandLineArguments& cmd);
    void FindNonPBCTriangles(CommandLineArguments& cmd);
    void CutMesh(CommandLineArguments& cmd);
    void ConvertToNonPBCMesh(CommandLineArguments& cmd);
    void ScaleMesh(CommandLineArguments& cmd);
};