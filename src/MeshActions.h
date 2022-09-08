#pragma once

#include <vector>
#include <string>
#include <array>
#include <memory>

#include "tools/Assert.h"
#include "MeshPlaneIntersection.hpp"
#include "tools/CommandLineArguments.h"
#include "Mesh.h"
#include "tools/CommonTypes.h"
#include "Curvature.h"
#include "tools/InputParser.h"
#include "MeshRefineStrategy.h"
#include "CurvatureCurveFit.h"
#include "CurvatureJetFit.h"
#include "tools/InputParser.h"
#include "tools/CommonOperations.h"
#include "tools/Algorithm.h"
#include "LinAlgTools.h"

namespace MeshActions
{
    using INT3 = std::array<int,3>;
    using INT2 = std::array<int,2>;
    using Real = CommonTypes::Real;
    using Real3= CommonTypes::Real3;
    using Real2= CommonTypes::Real2;
    using curveptr = std::unique_ptr<Curvature>;
    using refineptr= std::unique_ptr<MeshRefineStrategy>;

    // translate a mesh
    void TranslateMesh(CommandLineArguments& cmd);

    // calculate the mesh of a mesh using the method of curve fitting (second fundamental form)
    void CurveFit(CommandLineArguments& cmd);

    // calculate the curvature of a mesh using the method of jet fitting
    void JetFit(CommandLineArguments& cmd);

    // calculate the curvature of a mesh using the method of tensor
    void TensorFit(CommandLineArguments& cmd);

    // calculate the curvature of a mesh using the method of FDM
    void FDMFit(CommandLineArguments& cmd);

    // curvature evolution 
    void CurveEvolution(CommandLineArguments& cmd);

    // smooth a mesh using curvature flow (implicit fairing)
    void CurvatureFlow(CommandLineArguments& cmd);

    // function that finds the indices of non periodic triangles in a periodic mesh 
    void FindNonPBCTriangles(CommandLineArguments& cmd);

    // cut the mesh 
    void CutMesh(CommandLineArguments& cmd);

    // convert a pbc mesh to non pbc mesh
    void ConvertToNonPBCMesh(CommandLineArguments& cmd);

    // scale mesh by some number 
    void ScaleMesh(CommandLineArguments& cmd);

    // put color onto vertex based on curvature value
    void ColorVertex(CommandLineArguments& cmd);

    // project 3d curvature onto 2d using Moller Trumbore Ray-Triangle intersection algorithm
    void Project3dCurvature(CommandLineArguments& cmd);

    // function that constructs a new mesh such that none of the vertices is within a certain distance of the reference mesh 
    void MeshDistanceCutoff(CommandLineArguments& cmd);

    // Find the minimum distance of a mesh against another reference mesh
    void MinimumMeshDistance(CommandLineArguments& cmd);

    // find mesh plane intersection
    void MeshPlaneIntersection(CommandLineArguments& cmd);

    // find vertex neighbors 
    void FindVertexNeighbors(CommandLineArguments& cmd);

    // find boundary vertices 
    void FindBoundaryVertices(CommandLineArguments& cmd);

    // find the 2d projection of a 3d mesh viewed from a certain angle 
    void Project3dMesh(CommandLineArguments& cmd);

    // cut overlapped regions of 2 proteins based on Ray triangle intersection
    void CutOverlappedRegion(CommandLineArguments& cmd);
};