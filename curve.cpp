#include "src/Curvature.h"
#include "src/Mesh.h"
#include "tools/InputParser.h"

#include <string>
#include <vector>
#include <memory>

int main(int argc, char** argv)
{
    using curveptr = std::unique_ptr<Curvature>;
    using Meshptr  = std::unique_ptr<Mesh>;
    std::string fname = argv[1];

    InputParser ip;

    ParameterPack pack;
    ip.ParseFile(fname, pack);

    auto vertexPack = pack.findParamPack("vertexfile", ParameterPack::KeyType::Required);
    auto TrianglePack = pack.findParamPack("trianglefile", ParameterPack::KeyType::Required);
    auto curvaturePack = pack.findParamPack("curvature", ParameterPack::KeyType::Required);

    std::string vertexFileName;
    std::string triangleFileName;
    std::string curvaturetype;
    vertexPack->ReadString("name", ParameterPack::KeyType::Required, vertexFileName);
    TrianglePack->ReadString("name", ParameterPack::KeyType::Required, triangleFileName);
    curvaturePack->ReadString("type", ParameterPack::KeyType::Required,curvaturetype);


    Meshptr mesh_ = Meshptr(new Mesh(pack));    
    
    MeshTools::readPLY(vertexFileName, *mesh_);
    MeshTools::readPLYTriangle(triangleFileName, *mesh_);

    mesh_ -> refine();

    auto& vertices = mesh_->accessvertices();
    auto& triangles= mesh_->accesstriangles();

    // for (int i=0;i<mesh_.getNumTriangles();i++)
    // {
    //     auto& t = triangles[i];

    //     for (int j=0;j<3;j++)
    //     {
    //         std::cout << t.triangleindices_[j] << " ";
    //     }
    //     std::cout << "\n";
    // }

    // for (int i=0;i<mesh_.getNumVertices();i++)
    // {
    //     auto& v = vertices[i];

    //     for (int j=0;j<3;j++)
    //     {
    //         std::cout << v.position_[j] << " ";
    //     }

    //     std::cout << "\n";

    //     for (int j=0;j<3;j++)
    //     {
    //         std::cout << v.normals_[j] << " ";
    //     }
    //     std::cout << "\n";
    // }

    CurvatureInput input = { const_cast<ParameterPack&>(*curvaturePack), *mesh_};
    curveptr curve = curveptr(CurvatureRegistry::Factory::instance().create(curvaturetype, input));

    curve -> calculate();
    curve -> printOutput();
}