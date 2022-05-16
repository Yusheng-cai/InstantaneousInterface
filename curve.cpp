#include "src/Curvature.h"
#include "src/Mesh.h"
#include "tools/InputParser.h"
#include "src/Field.h"
#include "tools/CommonTypes.h"
#include "src/marching_cubes.hpp"
#include "src/Registry.h"
#include "src/MeshRefineStrategy.h"

#include <string>
#include <vector>
#include <memory>

int main(int argc, char** argv)
{
    using curveptr = std::unique_ptr<Curvature>;
    using refineptr= std::unique_ptr<MeshRefineStrategy>;
    using Meshptr  = std::unique_ptr<Mesh>;
    using fieldptr = std::unique_ptr<Field>;
    using Real3    = CommonTypes::Real3;
    using Real     = CommonTypes::Real;
    using Range    = std::array<Real,2>;
    std::string fname = argv[1];

    InputParser ip;
    std::vector<curveptr> curves_;
    std::vector<refineptr> refines_;

    ParameterPack pack;
    ip.ParseFile(fname, pack);

    Registry reg_;

    auto plyPack    = pack.findParamPack("plyfile", ParameterPack::KeyType::Optional);
    auto curvaturePack = pack.findParamPacks("curvature", ParameterPack::KeyType::Optional);
    auto refinePack = pack.findParamPacks("refine", ParameterPack::KeyType::Optional);

    // initialize the curvatures 
    for (int i=0;i<curvaturePack.size();i++)
    {
        auto& curvepack = curvaturePack[i];
        CurvatureInput input = { const_cast<ParameterPack&>(*curvepack)};
        std::string curvaturetype;

        curvepack -> ReadString("type", ParameterPack::KeyType::Required, curvaturetype);
        curveptr curve = curveptr(CurvatureRegistry::Factory::instance().create(curvaturetype, input));

        curves_.push_back(std::move(curve));

        reg_.registerCurvature(curves_[i]->getName(), *curves_[i]);
    }

    // initialize refinements
    for (int i=0;i<refinePack.size();i++)
    {
        auto& refinep = refinePack[i];
        MeshRefineStrategyInput input = {const_cast<ParameterPack&>(*refinep)};
        std::string refineType;

        refinep -> ReadString("type", ParameterPack::KeyType::Required, refineType);
        refineptr r = refineptr(MeshRefineStrategyFactory::Factory::instance().create(refineType, input));

        refines_.push_back(std::move(r));
        reg_.registerMeshRefinement(refines_[i]->getName(), *refines_[i]);
    }

    if (plyPack == nullptr)
    {
        std::cout << "No ply file is provided." << std::endl;
    }
    else
    {
        std::string plyFileName;
        plyPack->ReadString("name", ParameterPack::KeyType::Required,plyFileName);

        const auto meshPacks = pack.findParamPacks("Mesh", ParameterPack::KeyType::Optional);

        // special case when no mesh is provided
        if (meshPacks.size() == 0)  
        {
            Meshptr mesh_ = Meshptr(new Mesh(nullptr));    
        
            MeshTools::readPLYlibr(plyFileName, *mesh_);

            for (int i=0;i<refines_.size();i++)
            {
                refines_[i] -> refine(*mesh_);
            }

            mesh_ -> print();

            auto& vertices = mesh_->accessvertices();
            auto& triangles= mesh_->accesstriangles();


            for (int i=0;i<curvaturePack.size();i++)
            {
                auto& curve = curves_[i];
                curve -> calculate(*mesh_);
                curve -> printOutput();
            }
        }

        for (int i=0;i<meshPacks.size();i++)
        {
            std::cout << "Doing mesh " << i << std::endl;
            Meshptr mesh_ = Meshptr(new Mesh(meshPacks[i]));    
            std::vector<std::string> refinementNames;
            meshPacks[i] -> ReadVectorString("refinement", ParameterPack::KeyType::Optional, refinementNames);
            MeshTools::readPLYlibr(plyFileName, *mesh_);

            for (int i=0;i<refinementNames.size();i++)
            {
                std::string n = refinementNames[i];
                auto& refineMent = reg_.getMeshRefineStrat(n);

                refineMent.refine(*mesh_);
            }

            mesh_ -> print();

            auto& vertices = mesh_->accessvertices();
            auto& triangles= mesh_->accesstriangles();


            for (int i=0;i<curvaturePack.size();i++)
            {
                auto& curve = curves_[i];
                curve -> calculate(*mesh_);
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
        bool pbc=false;
        Range xrange, yrange, zrange;
        Real val;
        std::string fname;
        MarchingCubes mwrapper_;
        fpack -> ReadArrayNumber("xrange", ParameterPack::KeyType::Required, xrange);
        fpack -> ReadArrayNumber("yrange", ParameterPack::KeyType::Required, yrange);
        fpack -> ReadArrayNumber("zrange", ParameterPack::KeyType::Required, zrange);
        fpack -> ReadArrayNumber("dimension", ParameterPack::KeyType::Required, index);
        fpack -> ReadNumber("file", ParameterPack::KeyType::Required,fname);
        fpack -> ReadNumber("isosurfacevalue", ParameterPack::KeyType::Required,val);
        fpack -> Readbool("pbc", ParameterPack::KeyType::Optional, pbc);

        // find the field pointer
        fieldptr f = fieldptr(new Field(index[0], index[1],index[2], xrange, yrange, zrange));
        FieldTools::FieldReader(fname, *f);

        auto mpack = fpack -> findParamPack("Mesh", ParameterPack::KeyType::Required);
        Meshptr mesh_ = Meshptr(new Mesh(mpack));
        std::vector<std::string> refinements;
        if (mpack != nullptr)
        {
            mpack -> ReadVectorString("refinement", ParameterPack::KeyType::Optional, refinements);
        }

        // calculate the marching cubes
        mwrapper_.triangulate_field(*f, *mesh_, val, pbc);
        auto cpack = fpack -> findParamPacks("curvature", ParameterPack::KeyType::Optional);

        for (int i=0;i<refinements.size();i++)
        {
            std::string name = refinements[i];
            auto& r = reg_.getMeshRefineStrat(name);
            r.refine(*mesh_);
        }

        mesh_ -> print();

        for (int j=0;j<cpack.size();j++)
        {
            std::string curvetype;
            auto cp = cpack[j];
            cp -> ReadString("type", ParameterPack::KeyType::Required, curvetype);

            CurvatureInput input = { const_cast<ParameterPack&>(*cp)};
            curveptr c = curveptr(CurvatureRegistry::Factory::instance().create(curvetype, input));

            c -> calculate(*mesh_);
            c -> printOutput();
        }
    }
}