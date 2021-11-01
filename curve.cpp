#include "src/Curvature.h"
#include "src/Mesh.h"
#include "tools/InputParser.h"
#include "src/Field.h"
#include "tools/CommonTypes.h"
#include "src/MarchingCubesWrapper.h"

#include <string>
#include <vector>
#include <memory>

int main(int argc, char** argv)
{
    using curveptr = std::unique_ptr<Curvature>;
    using Meshptr  = std::unique_ptr<Mesh>;
    using fieldptr = std::unique_ptr<Field>;
    using Real3    = CommonTypes::Real3;
    using Real     = CommonTypes::Real;
    using Range    = std::array<Real,2>;
    std::string fname = argv[1];

    InputParser ip;

    ParameterPack pack;
    ip.ParseFile(fname, pack);

    auto plyPack    = pack.findParamPack("plyfile", ParameterPack::KeyType::Optional);
    auto curvaturePack = pack.findParamPacks("curvature", ParameterPack::KeyType::Optional);

    if (plyPack == nullptr)
    {
        std::cout << "No ply file is provided." << std::endl;
    }
    else
    {
        std::string vertexFileName;
        std::string triangleFileName;
        std::string plyFileName;
        plyPack->ReadString("name", ParameterPack::KeyType::Required,plyFileName);

        const auto meshPacks = pack.findParamPacks("Mesh", ParameterPack::KeyType::Optional);

        // special case when no mesh is provided
        if (meshPacks.size() == 0)  
        {
            Meshptr mesh_ = Meshptr(new Mesh(nullptr));    
        
            MeshTools::readPLYlibr(plyFileName, *mesh_);

            mesh_ -> refine();
            mesh_ -> CalcTriangleAreaAndFacetNormals();
            mesh_ -> print();

            auto& vertices = mesh_->accessvertices();
            auto& triangles= mesh_->accesstriangles();


            for (int i=0;i<curvaturePack.size();i++)
            {
                auto& curvepack = curvaturePack[i];
                CurvatureInput input = { const_cast<ParameterPack&>(*curvepack), *mesh_};
                std::string curvaturetype;

                curvepack -> ReadString("type", ParameterPack::KeyType::Required, curvaturetype);
                curveptr curve = curveptr(CurvatureRegistry::Factory::instance().create(curvaturetype, input));

                curve -> calculate();
                curve -> printOutput();
            }

        }

        for (int i=0;i<meshPacks.size();i++)
        {
            std::cout << "Doing mesh " << i << std::endl;
            Meshptr mesh_ = Meshptr(new Mesh(meshPacks[i]));    
        
            MeshTools::readPLYlibr(plyFileName, *mesh_);

            mesh_ -> refine();
            mesh_ -> CalcTriangleAreaAndFacetNormals();
            mesh_ -> print();

            auto& vertices = mesh_->accessvertices();
            auto& triangles= mesh_->accesstriangles();


            for (int i=0;i<curvaturePack.size();i++)
            {
                auto& curvepack = curvaturePack[i];
                CurvatureInput input = { const_cast<ParameterPack&>(*curvepack), *mesh_};
                std::string curvaturetype;

                curvepack -> ReadString("type", ParameterPack::KeyType::Required, curvaturetype);
                curveptr curve = curveptr(CurvatureRegistry::Factory::instance().create(curvaturetype, input));

                curve -> calculate();
                curve -> printOutput();
            }
        }
    }

    // find the field
    auto fieldPack = pack.findParamPacks("field", ParameterPack::KeyType::Optional);

    if (fieldPack.size() == 0)
    {
        std::cout << "No field is provided, program will now exit .... " << std::endl;
    }

    for (int i=0;i<fieldPack.size();i++)
    {
        auto fpack = fieldPack[i];

        std::array<int,3> index;
        int Nx, Ny, Nz;
        Range xrange, yrange, zrange;
        Real val;
        std::string fname;
        MarchingCubesWrapper mwrapper_;
        fpack -> ReadArrayNumber("xrange", ParameterPack::KeyType::Required, xrange);
        fpack -> ReadArrayNumber("yrange", ParameterPack::KeyType::Required, yrange);
        fpack -> ReadArrayNumber("zrange", ParameterPack::KeyType::Required, zrange);
        fpack -> ReadArrayNumber("dimension", ParameterPack::KeyType::Required, index);
        fpack -> ReadNumber("file", ParameterPack::KeyType::Required,fname);
        fpack -> ReadNumber("isosurfaceval", ParameterPack::KeyType::Required,val);

        // find the field pointer
        fieldptr f = fieldptr(new Field(index[0], index[1],index[2], xrange, yrange, zrange));
        FieldTools::FieldReader(fname, *f);

        auto mpack = fpack -> findParamPack("Mesh", ParameterPack::KeyType::Required);
        Meshptr mesh_ = Meshptr(new Mesh(mpack));

        // calculate the marching cubes
        mwrapper_.calculate(*f, *mesh_, val);
        auto cpack = fpack -> findParamPacks("curvature", ParameterPack::KeyType::Optional);

        mesh_ -> refine();
        mesh_ -> print();

        for (int j=0;j<cpack.size();j++)
        {
            std::string curvetype;
            auto cp = cpack[j];
            cp -> ReadString("type", ParameterPack::KeyType::Required, curvetype);

            CurvatureInput input = { const_cast<ParameterPack&>(*cp), *mesh_};
            curveptr c = curveptr(CurvatureRegistry::Factory::instance().create(curvetype, input));

            c -> calculate();
            c -> printOutput();
        }
    }
}