#include "AverageField.h"

namespace DensityFieldRegistry
{
    registry_<AverageField> registerAverageField("averagefield");
}

AverageField::AverageField(const DensityFieldInput& input)
:DensityField(input)
{}

void AverageField::calculate()
{
    auto atomgroup = getAtomGroup(atomGroupName_);
    auto& atoms = atomgroup.getAtoms();
    field_.zero();

    for (int i=0;i<atoms.size();i++)
    {
        if( bound_box_ -> isInside(atoms[i].position))
        {
            std::cout << "Position before correction = " << atoms[i].position[0] << " " << atoms[i].position[1] << " " << \
            atoms[i].position[2] << std::endl;

            Real3 correctedPos = bound_box_->PutInBoundingBox(atoms[i].position);

            std::cout << "Positon after correction = " << correctedPos[0] << " " << correctedPos[1] << " " << correctedPos[2] << std::endl;
            index3 Index       = field_.getClosestGridIndex(correctedPos);
            std::cout << "Closest position = " << Index[0] << " " << Index[1] << " " << Index[2] << std::endl;

            for(int j=0;j<offsetIndex_.size();j++)
            {
                index3 RealIndex;
                for (int k=0;k<3;k++)
                {
                    RealIndex[k] = Index[k] + offsetIndex_[j][k];
                }

                std::cout << "Index = " << RealIndex[0] << " " << RealIndex[1] << " " << RealIndex[2] << std::endl;
                Real3 latticepos = field_.getPositionOnGrid(RealIndex[0], RealIndex[1], RealIndex[2]);
                std::cout << "Position on grid = " << latticepos[0] << " " << latticepos[1] << " " << latticepos[2] << std::endl;

                Real3 diff;
                diff.fill(0);

                for (int k=0;k<3;k++)
                {
                    diff[k] = std::abs(latticepos[k] - correctedPos[k]);
                }

                Real val = GaussianCoarseGrain::GaussianCoarseGrainFunction(diff, sigma_);

                std::cout << "The value of the field = " << val << std::endl;

                field_(RealIndex[0], RealIndex[1], RealIndex[2]) += val;
            }
        }
    }    
}