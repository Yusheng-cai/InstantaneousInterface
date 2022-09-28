#include "MeshGeneration.h"

MeshGeneration::MeshGeneration(ParameterPack& pack)
: pack_(pack)
{
    pack_.ReadVectorString("outputs", ParameterPack::KeyType::Optional, OutputNames_);
    pack_.ReadVectorString("outputNames", ParameterPack::KeyType::Optional, OutputFileNames_);

    output_.registerOutputFunc("ply", [this](std::string name)-> void {this -> printPLY(name);});
}

void MeshGeneration::printPLY(std::string name)
{
    MeshTools::writePLY(name, mesh_);
}

void MeshGeneration::printOutputs()
{
    for (int i=0;i<OutputNames_.size();i++){
        output_.getOutputFuncByName(OutputNames_[i])(OutputFileNames_[i]);
    }
}