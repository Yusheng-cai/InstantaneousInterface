#include "Curvature.h"

Curvature::Curvature(CurvatureInput& input)
:mesh_(input.mesh_)
{
    bool read = input.pack.ReadString("output", ParameterPack::KeyType::Optional, OutputName_);

    if (read)
    {
        ofs_.open(OutputName_);
        ASSERT((ofs_.is_open()), "The file with name " << OutputName_ << " is not opened.");
    }
}

void Curvature::close()
{
    if (ofs_.is_open())
    {
        ofs_.close();
    }
}