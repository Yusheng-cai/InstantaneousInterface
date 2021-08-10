#include "AverageField.h"

namespace DensityFieldRegistry
{
    registry_<AverageField> registerAverageField("averagefield");
}

AverageField::AverageField(const DensityFieldInput& input)
:DensityField(input)
{}

void AverageField::precalculateDensity()
{
    DensityPreCalculate_.clear(); 
    DensityPreCalculate_.resize(offsetIndex_.size());

    Real3 originPos = field_.getPositionOnGrid(0,0,0);

    #pragma omp parallel for
    for (int i=0;i<offsetIndex_.size();i++)
    {
        index3 index = offsetIndex_[i];

        Real3 pos = field_.getPositionOnGrid(index[0], index[1], index[2]);

        Real3 distance;
        Real sq_dist;
        simstate_.getSimulationBox().calculateDistance(pos, originPos, distance, sq_dist);

        Real val = GaussianCoarseGrain::GaussianCoarseGrainFunction(distance, sigma_);
        DensityPreCalculate_[i] = val;
    }
}

void AverageField::calculate()
{
    precalculateDensity();

    auto atomgroup = getAtomGroup(atomGroupName_);
    auto& atoms = atomgroup.getAtoms();

    // set up the omp buffers
    FieldBuffer_.set_master_object(field_);
    AtomIndicesInside_.clear();
    AtomIndicesBuffer_.set_master_object(AtomIndicesInside_);

    #pragma omp parallel
    {
        auto& indices_buffer = AtomIndicesBuffer_.access_buffer_by_id();
        #pragma omp for
        for(int i=0;i<atoms.size();i++)
        {
            if (bound_box_ -> isInside(atoms[i].position))
            {
                indices_buffer.push_back(i);
            }
        }
    }

    int size = AtomIndicesInside_.size();
    for (auto it = AtomIndicesBuffer_.beginworker();it != AtomIndicesBuffer_.endworker();it++)
    {
        size += it -> size();
    }

    AtomIndicesInside_.reserve(size);

    for (auto it = AtomIndicesBuffer_.beginworker();it != AtomIndicesBuffer_.endworker();it++)
    {
        AtomIndicesInside_.insert(AtomIndicesInside_.end(), it->begin(), it -> end());
    }

    #pragma omp parallel
    {
        auto& fieldbuf = FieldBuffer_.access_buffer_by_id();

        #pragma omp for simd
        for (int i=0;i<AtomIndicesInside_.size();i++)
        {
            int indices = AtomIndicesInside_[i];
            //std::cout << "Position before correction = " << atoms[indices].position[0] << " " << atoms[indices].position[1] << " " << \
            atoms[indices].position[2] << std::endl;

            Real3 correctedPos = bound_box_->PutInBoundingBox(atoms[indices].position);

            //std::cout << "Positon after correction = " << correctedPos[0] << " " << correctedPos[1] << " " << correctedPos[2] << std::endl;
            index3 Index       = fieldbuf.getClosestGridIndex(correctedPos);
            //std::cout << "Closest position = " << Index[0] << " " << Index[1] << " " << Index[2] << std::endl;

            for(int j=0;j<offsetIndex_.size();j++)
            {
                index3 RealIndex;
                for (int k=0;k<3;k++)
                {
                    RealIndex[k] = Index[k] + offsetIndex_[j][k];
                }

                // std::cout << "Index = " << RealIndex[0] << " " << RealIndex[1] << " " << RealIndex[2] << std::endl;
                Real3 latticepos = fieldbuf.getPositionOnGrid(RealIndex[0], RealIndex[1], RealIndex[2]);
                // std::cout << "Position on grid = " << latticepos[0] << " " << latticepos[1] << " " << latticepos[2] << std::endl;

                // Real3 diff;
                // diff.fill(0);

                // for (int k=0;k<3;k++)
                // {
                //     diff[k] = std::abs(latticepos[k] - correctedPos[k]);
                // }

                // Real val = GaussianCoarseGrain::GaussianCoarseGrainFunction(diff, sigma_);

                Real val = DensityPreCalculate_[j];
                // std::cout << "val = " << val << std::endl;
                // std::cout << "index = " << RealIndex[0] << " " << RealIndex[1] << " " << RealIndex[2] << std::endl;

                // std::cout << "The value of the field = " << val << std::endl;

                fieldbuf(RealIndex[0], RealIndex[1], RealIndex[2]) += val;
            }
        }    
    }

    #pragma omp parallel for
    for (int i=0;i<field_.totalSize();i++)
    {
        for (auto it = FieldBuffer_.beginworker();it != FieldBuffer_.endworker();it++)
        {
            field_.accessField()[i] += it->accessField()[i];
        }
    }
}

void AverageField::finishCalculate()
{
    MarchingCubes_.calculate(field_, vertices_, triangles_, isoSurfaceVal_);

    // for (int i=0;i<vertices_.size();i++)
    // {
    //     for (int j=0;j<3;j++)
    //     {
    //         std::cout << vertices_[i].position_[j] << " ";
    //     }
    //     std::cout << "\n";
    // }
}

void AverageField::printOutputIfOnStep()
{
    index3 N = field_.getN();
    int Nx = N[0];
    int Ny = N[1];
    int Nz = N[2];

    for (int i=0;i<field_.accessField().size();i++) 
    {
        ofs_ << field_.accessField()[i];

        ofs_ << " ";
    }
}