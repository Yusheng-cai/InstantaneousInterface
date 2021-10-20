#include "Curvature.h"

Curvature::Curvature(CurvatureInput& input)
:mesh_(input.mesh_)
{
    input.pack.ReadVectorString("outputs", ParameterPack::KeyType::Optional, OutputNames_);
    input.pack.ReadVectorString("outputNames", ParameterPack::KeyType::Optional, OutputFileNames_);
    outputs_.registerOutputFunc("curvature", [this](std::string name) -> void { this -> printCurvature(name);});

    ASSERT(( OutputNames_.size() == OutputFileNames_.size()), "The number of outputs does not agree with number of output files.");
}

void Curvature::printOutput()
{
    for (int i=0;i<OutputNames_.size();i++)
    {
        outputs_.getOutputFuncByName(OutputNames_[i])(OutputFileNames_[i]);
    }
}