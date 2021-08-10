#include "AverageField.h"

namespace DensityFieldRegistry
{
    registry_<AverageField> registerAverageField("averagefield");
}

AverageField::AverageField(const DensityFieldInput& input)
:DensityField(input)
{
    // Read in the output file
    bool outputread = input.pack_.ReadString("field_output", ParameterPack::KeyType::Optional, fieldOutputFileName_);
    if (outputread)
    {
        fieldofs_.open(fieldOutputFileName_);
        ASSERT((fieldofs_.is_open()), "The file with name " << fieldOutputFileName_ << " is not opened.");
    }

    bool triangleoutput = input.pack_.ReadString("triangle_output", ParameterPack::KeyType::Optional, triangleIndicesFileName_);
    if (triangleoutput)
    {
        triangleIndicesofs_.open(triangleIndicesFileName_);
        ASSERT((triangleIndicesofs_.is_open()), "The file with name " << triangleIndicesFileName_ << " is not opened.");
    }

    bool vertexoutput = input.pack_.ReadString("vertex_output", ParameterPack::KeyType::Optional, vertexFileName_);
    if (vertexoutput)
    {
        vertexofs_.open(vertexFileName_);
        ASSERT((vertexofs_.is_open()), "The file with name " << vertexFileName_ << " is not opened.");
    }

    bool vertexNormalOutput = input.pack_.ReadString("normal_output", ParameterPack::KeyType::Optional, vertexNormalFileName_);
    if (vertexNormalOutput)
    {
        vertexNormalofs_.open(vertexNormalFileName_);
        ASSERT((vertexNormalofs_.is_open()), "The file with name " << vertexNormalFileName_ << " is not opened.");
    }
}

void AverageField::calculate()
{
    if (PreCalculateDensity_)
    {
        CalculateUsingPreCalculatedDensity(); 
    }
    else
    {
        CalculateOnTheFly();
    }
}

void AverageField::finishCalculate()
{
    auto& fieldVec = field_.accessField();
    std::cout << "Total frames to be calculated = " << simstate_.getTotalFramesToBeCalculated() << std::endl;
    Real avgFac = 1.0/simstate_.getTotalFramesToBeCalculated();

    // performing average on the field
    for (int i=0;i<fieldVec.size();i++)
    {
        fieldVec[i] = avgFac*fieldVec[i]; 
    }

    MarchingCubes_.calculate(field_, vertices_, triangles_, isoSurfaceVal_);
}

void AverageField::printOutputIfOnStep()
{}

void AverageField::printFinalOutput()
{
    index3 N = field_.getN();
    int Nx = N[0];
    int Ny = N[1];
    int Nz = N[2];

    if (fieldofs_.is_open())
    {
        for (int i=0;i<field_.accessField().size();i++) 
        {
            fieldofs_ << field_.accessField()[i];

            fieldofs_ << " ";
        }
        fieldofs_.close();
    }


    if (vertexofs_.is_open())
    {
        for (int i=0;i<vertices_.size();i++)
        {
            for (int j=0;j<3;j++)
            {
                vertexofs_ << vertices_[i].position_[j];

                vertexofs_ << " ";
            }

            vertexofs_ << "\n";
        }
        vertexofs_.close();
    }

    if (vertexNormalofs_.is_open())
    {
        for (int i=0;i<vertices_.size();i++)
        {
            for (int j=0;j<3;j++)
            {
                vertexNormalofs_ << vertices_[i].normals_[j];
                
                vertexNormalofs_ << " ";
            }

            vertexNormalofs_ << "\n";
        }
        vertexNormalofs_.close();
    }

    if (triangleIndicesofs_.is_open())
    {
        for (int i=0;i<triangles_.size();i++)
        {
            for (int j=0;j<3;j++)
            {
                triangleIndicesofs_ << triangles_[i].triangleindices_[j];

                triangleIndicesofs_ << " ";
            }
            triangleIndicesofs_ << "\n";
        }

        triangleIndicesofs_.close();
    }

}