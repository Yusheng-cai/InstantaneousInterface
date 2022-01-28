#include "src/MeshGeneration.h"
#include "tools/Assert.h"
#include "tools/InputParser.h"

#include <vector>
#include <string>
#include <array>
#include <memory>

int main(int argc, char** argv)
{
    using meshgen = std::unique_ptr<MeshGeneration>; 

    ASSERT((argc==2), "The input file is not provided.");
    std::string fname = argv[1];

    InputParser parser;
    ParameterPack pack;
    parser.ParseFile(fname, pack);

    const auto meshGenPacks = pack.findParamPacks("meshgen", ParameterPack::KeyType::Optional);

    for (auto m : meshGenPacks)
    {
        std::string type;
        m->ReadString("type", ParameterPack::KeyType::Required, type);
        meshgen gen = meshgen(MeshGenerationFactory::factory::instance().create(type, const_cast<ParameterPack&>(*m)));

        gen->generate();

        gen->printOutputs();
    }
};