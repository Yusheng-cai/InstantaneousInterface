#include <iostream>
#include <vector>
#include <array>
#include <string>
#include <map>
#include <functional>

#include "tools/CommandLineArguments.h"
#include "tools/CommonTypes.h"
#include "tools/Assert.h"
#include "src/Curvature.h"
#include "src/Mesh.h"
#include "src/MeshActions.h"

using ActionFunction = std::function<void(CommandLineArguments& cmd)>;
using mapFunction    = std::map<std::string, ActionFunction>;
using mapUsage       = std::map<std::string, std::string>;

void RegisterAction(std::string name, std::string usage, ActionFunction func, mapFunction& mapF, mapUsage& mapU);
void RegisterAllActions(mapFunction& mapF, mapUsage& mapU);

int main(int argc, char** argv){
    CommandLineArguments cmd(argc, argv);

    std::string help_key_ = "help";
    std::string operation_;
    mapFunction MapNameToAction;
    mapUsage    MapNameToUsage;

    RegisterAllActions(MapNameToAction, MapNameToUsage);

    if (cmd.has_key("help") && ! cmd.has_key("op")){
        std::cout << "Usage: mesh_op -op OPERATION " << "\n";
        std::cout << "operation  Usage" << "\n";
        for (auto it = MapNameToUsage.begin(); it != MapNameToUsage.end(); it ++){
            std::cout << it -> first << "\t" << it -> second << "\n";
        }
        return 0;
    }


    // find the operation
    cmd.readString("op",CommandLineArguments::Keys::Required,operation_);

    if (cmd.has_key("help") && cmd.has_key("op")){
        auto it = MapNameToUsage.find(operation_);
        std::cout << it -> first << "\t" << it -> second << "\n";

        return 0;
    }

    // perform the operation
    auto it = MapNameToAction.find(operation_);
    ASSERT((it != MapNameToAction.end()), "Operation " << operation_ << " is not registered.");

    it -> second(cmd);

    return 0;
}

void RegisterAction(std::string name, std::string usage, ActionFunction func, mapFunction& mapF, mapUsage& mapU)
{
    auto Fit = mapF.find(name);
    ASSERT((Fit == mapF.end()), "The operation " << name << " is registered more than once.");

    auto Uit = mapU.find(name);
    ASSERT((Uit == mapU.end()), "The operation " << name << " is registered more than once.");

    mapF.insert(std::make_pair(name,func));
    mapU.insert(std::make_pair(name, usage));
}

void RegisterAllActions(mapFunction& mapF, mapUsage& mapU)
{
    RegisterAction("translate", "-op translate -trans dx dy dz -i input.ply -o output.ply", \
                    [](CommandLineArguments& cmd)-> void {return MeshActions::TranslateMesh(cmd);}, \
                    mapF, mapU);
    RegisterAction("curvefit", "-op curvefit -box[optional] x y z -i input.ply -o curve.out -neighbors number[int]", \
                    [](CommandLineArguments& cmd)-> void {return MeshActions::CurveFit(cmd);}, \
                    mapF, mapU);
    RegisterAction("jetfit", "-op jetfit -box[optional] x y z -i input.ply -o curve.out -fc faceC.out -neighbors number[int] -MongeCoefficient number[int] -degree number[int]",\
                              [](CommandLineArguments& cmd)->void{return MeshActions::JetFit(cmd);}, \
                              mapF, mapU);
    RegisterAction("quadraticfit", "-op quadraticfit -box[optional] x y z -i input.ply -o quadraticfit.out -neighbors number[int] -MonteCarlo bool -MonteCarloN number", \
                              [](CommandLineArguments& cmd)->void{return MeshActions::QuadraticCurveFit(cmd);}, \
                              mapF, mapU);
    RegisterAction("tensorfit", "-op tensorfit -box[optional] x y z -i input.ply -o tensor.out -fc tensor_faceC.out", \
                              [](CommandLineArguments& cmd)-> void {return MeshActions::TensorFit(cmd);}, \
                              mapF, mapU);
    RegisterAction("FDMfit", "-op FDMfit -box[optional] x y z -i input.ply -o FDMfit.out -fc FDMfit_faceC.out", \
                              [](CommandLineArguments& cmd)->void {return MeshActions::FDMFit(cmd);}, \
                              mapF, mapU);
    RegisterAction("curvatureflow", "-op curvatureflow -box[optional] x y z -i input.ply -o output.ply -pbcOutput[optional] true/false -Decimate[optional] true -iteration number[int] -lambdadt number[float]", \
                    [](CommandLineArguments& cmd)-> void {return MeshActions::CurvatureFlow(cmd);}, \
                    mapF, mapU);
    RegisterAction("NonPBCFace", "-op NonPBCFace -box[required] x y z -i input.ply -o NonPBC.out",\
                    [](CommandLineArguments& cmd) -> void {return MeshActions::FindNonPBCTriangles(cmd);}, \
                    mapF, mapU);
    RegisterAction("cut", "-op cut -volume x y z -i input.ply -o output.ply", \
                    [](CommandLineArguments& cmd)-> void {return MeshActions::CutMesh(cmd);}, \
                    mapF, mapU);
    RegisterAction("ConvertToNonPBCMesh", "-op ConvertToNonPBCMesh -box[required] x y z -i[required] input.ply -o [nonpbc.ply]", \
                    [](CommandLineArguments& cmd) -> void {return MeshActions::ConvertToNonPBCMesh(cmd);}, \
                    mapF, mapU);
    RegisterAction("scale", "-op scale -scale[float] num -i input.ply -o scaled.ply", \
                    [](CommandLineArguments& cmd)-> void {return MeshActions::ScaleMesh(cmd);}, \
                    mapF, mapU);
    RegisterAction("ProjectCurvature", "-op ProjectCurvature -i input.ply -o curvature.out -ProjectedIndex[array2] dim1 dim2 -height[float] num -n n1 n2 -L L1 L2 -FaceCurvatureFile[required] fc.out -box[optional] x y z -RayDirection x y z", 
                                        [](CommandLineArguments& cmd)-> void {return MeshActions::Project3dCurvature(cmd);}, 
                                        mapF, mapU);
    RegisterAction("DistanceCutoff", "-op DistanceCutoff -i input.ply -ref ref.ply -o[optional] cut.ply -cutoff[float] distance \
                                      -box[optional] x y z",\
                                      [](CommandLineArguments& cmd)-> void {return MeshActions::MeshDistanceCutoff(cmd);}, \
                                      mapF, mapU);
    RegisterAction("ColorVertex", "-op ColorVertex -i input.ply -o output.ply -curvature[string] c.out -col[int] column -vmin[optional] \
                                   -vmax[optional]", \
                                    [](CommandLineArguments& cmd)-> void {return MeshActions::ColorVertex(cmd);}, \
                                    mapF, mapU);
    RegisterAction("MeshMinDistance", "-op MeshMinDistance -i input.ply -ref ref.ply -o output.out -box[optional] x y z", \
                                    [](CommandLineArguments& cmd)-> void {return MeshActions::MinimumMeshDistance(cmd);}, \
                                    mapF, mapU);
    RegisterAction("MeshPlaneIntersect", "-op MeshPlaneIntersect -i input.ply -o intersect.out -plane[Real3] x y z -point[Real3] x y z", \
                                    [](CommandLineArguments& cmd)-> void {return MeshActions::MeshPlaneIntersection(cmd);},\
                                    mapF, mapU);

    RegisterAction("BoundaryVertices", "-op BoundaryVertices -i input.ply -o boundary.out",\
                                    [](CommandLineArguments& cmd) -> void {return MeshActions::FindBoundaryVertices(cmd);},\
                                    mapF, mapU);
    RegisterAction("ProjectMesh", "-op ProjectMesh -i[vector] a.ply b.ply .. -L L1 L2 -n n1 n2 -RayDirection -ProjectedIndex ind1 ind2 -height h -box[optional] x y z",\
                                    [](CommandLineArguments& cmd)-> void {return MeshActions::Project3dMesh(cmd);}, \
                                    mapF, mapU);
    RegisterAction("CutOverlappedRegion", "-op CutOverlappedRegion -i input.ply -ref ref.ply -o output.ply -RayDirection -box[optional] x y z", \
                                           [](CommandLineArguments& cmd) -> void {return MeshActions::CutOverlappedRegion(cmd);}, \
                                           mapF, mapU);
    RegisterAction("DecimateDegenTriangle", "-op DecimateDegenTriangle -i input.ply -box[optional] x y z -o output.ply", \
                                           [](CommandLineArguments& cmd) -> void {return MeshActions::DecimateDegenerateTriangles(cmd);}, \
                                           mapF, mapU);
    RegisterAction("CurvatureEvolution", "-op CurvatureEvolution -i input.ply -box[optional] x y z -o output.ply -k0 kappa -stepsize size -maxiter max -neighbors[int] -tol tolerance", \
                                           [](CommandLineArguments& cmd) -> void {return MeshActions::CurvatureEvolution(cmd);}, \
                                           mapF, mapU);
    RegisterAction("FindIsolatedFace", "-op FindIsolatedFace -i input.ply -o isolated.out", \
                                           [](CommandLineArguments& cmd) -> void {return MeshActions::FindIsolatedFace(cmd);},\
                                           mapF, mapU);
    RegisterAction("CutTeethlikeFace", "-op CutTeethlikeFace -i input.ply -o cut.ply", \
                                           [](CommandLineArguments& cmd) -> void {return MeshActions::CutTeethlikeFace(cmd);},\
                                           mapF, mapU);
    RegisterAction("ConvertToSTL", "-op ConvertToSTL -i input.ply -o output.stl", \
                                           [](CommandLineArguments& cmd ) -> void {return MeshActions::ConvertToStL(cmd);}, \
                                           mapF, mapU);
    RegisterAction("ReplicatePBCMesh", "-op ReplicatePBCMesh -i input.ply -o output.ply -box x y z", \
                                           [](CommandLineArguments& cmd) -> void {return MeshActions::ReplicatePeriodicMesh(cmd);}, \
                                           mapF, mapU);
    RegisterAction("FindNeighborIndices", "-op FindNeighborIndices -i input.ply -n num -o output.out", \
                                           [](CommandLineArguments& cmd) -> void {return MeshActions::FindVertexNeighbors(cmd);}, \
                                           mapF, mapU);
    RegisterAction("AngleDistribution", "-op AngleDistribution -i input.ply -numbins num -box x y z -o dist.out", \
                                           [](CommandLineArguments& cmd) -> void {return MeshActions::TriangleAngleDistribution(cmd);}, \
                                           mapF, mapU);
    RegisterAction("SideLengthsDistribution", "-op SideLengthsDistribution -i input.ply -o sidelength.out -box x y z", \
                                           [](CommandLineArguments& cmd) -> void {return MeshActions::SideLengthsDistribution(cmd);}, \
                                           mapF, mapU);
    RegisterAction("MeshCleanup", "-op MeanCleanup -i input.ply -box x y z -o cleaned.ply -iterations num -edgeLengthCutoff num", \
                                           [](CommandLineArguments& cmd) -> void {return MeshActions::MeshCleanup(cmd);}, \
                                           mapF, mapU);
    RegisterAction("FlattenMesh", "-op FlattenMesh -i input.ply -o flat.ply -index num", \
                                           [](CommandLineArguments& cmd) -> void {return MeshActions::FlattenMesh(cmd);}, \
                                           mapF, mapU);
    RegisterAction("ConformingTriangles", "-op ConformingTriangles -box x y z", \
                                           [](CommandLineArguments& cmd) -> void {return MeshActions::ConformingTriangulations(cmd);}, \
                                           mapF, mapU);
    RegisterAction("SplitLongEdge", "-op SplitLongEdge -i input.ply -max_length m -o split.ply", \
                                           [](CommandLineArguments& cmd) -> void {return MeshActions::SplitLongEdges(cmd);}, \
                                           mapF, mapU);
    RegisterAction("ShiftMeshWithRef", "-op ShiftMeshWithRef -i input.ply -ref ref.ply -o out.ply", \
                                           [](CommandLineArguments& cmd) -> void {return MeshActions::ShiftMeshWithRef(cmd);}, \
                                           mapF, mapU);
    RegisterAction("ViewMeshWithData", "-op ViewMeshWithData -i input.ply -data data.out -col[int] -numSteps[optional] -min[optional] -max[optional]",\
                                           [](CommandLineArguments& cmd) -> void {return MeshActions::ViewMeshWithData(cmd);}, \
                                           mapF, mapU);
    RegisterAction("ViewMesh", "-op ViewMesh -i input.ply", \
                                           [](CommandLineArguments& cmd) -> void {return MeshActions::ViewMesh(cmd);}, \
                                           mapF, mapU);
    RegisterAction("FlattenMeshDimension", "-op FlattenMeshDimension -i input.ply -dim [int] -set [float] -o flat.ply" , \
                                           [](CommandLineArguments& cmd) -> void {return MeshActions::FlattenMeshDimension(cmd);}, \
                                           mapF, mapU);
    RegisterAction("DistanceBetweenMeshesMT", "-op DistanceBetweenMeshesMR -i input1.ply input2.ply -Ray [vec3] -box [vec3] -o output.out", \
                                           [](CommandLineArguments& cmd) -> void {return MeshActions::DistanceBetweenMeshesMT(cmd);}, \
                                           mapF, mapU);
    RegisterAction("IterativeClosestPoint", "-op IterativeClosestPoint -i input1.ply input2.ply -box [vec3] -o output.out", \
                                           [](CommandLineArguments& cmd) -> void {return MeshActions::IterativeClosestPoint(cmd);}, \
                                           mapF, mapU);
}
