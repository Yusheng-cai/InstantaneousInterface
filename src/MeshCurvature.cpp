#include "MeshCurvature.h"

namespace MeshRefineStrategyFactory
{
    registry_<MeshCurvature> registercurve("curvature");
}

MeshCurvature::MeshCurvature(MeshRefineStrategyInput& input)
: MeshRefineStrategy(input)
{
    // initialize curvature 
    auto cPack = pack_.findParamPack("Curvature", ParameterPack::KeyType::Required);
    cPack->ReadString("type", ParameterPack::KeyType::Required, curvaturetype_);
    CurvatureInput in = {const_cast<ParameterPack&>(*cPack)};
    curvatureCalc_ = curveptr(CurvatureRegistry::Factory::instance().create(curvaturetype_,in));

    // read in the necessary ones for update 
    pack_.ReadNumber("stepsize", ParameterPack::KeyType::Required, StepSize_);
    pack_.ReadNumber("tolerance", ParameterPack::KeyType::Optional, tol_);
    pack_.ReadNumber("k0", ParameterPack::KeyType::Required, meanCurvature_);
    pack_.ReadNumber("maxstep", ParameterPack::KeyType::Optional, maxStep);
    pack_.Readbool("FixBoundary", ParameterPack::KeyType::Optional,fixBoundary_);
}


void MeshCurvature::findVertices()
{
    const auto& vertices = mesh_->getvertices();

    if (fixBoundary_)
    {
        MeshTools::MapEdgeToFace(*mesh_, MapEdgeToFace_, MapVertexToEdge_);
        MeshTools::CalculateBoundaryVertices(*mesh_, MapEdgeToFace_, boundaryIndicator_);

        for (int i=0;i<vertices.size();i++)
        {
            if (! MeshTools::IsBoundary(i, boundaryIndicator_))
            {
                OutsideIndices_.push_back(i);
            }
        }

        // Find the unique indices 
        auto ip = std::unique(OutsideIndices_.begin(), OutsideIndices_.end());
        OutsideIndices_.resize(std::distance(OutsideIndices_.begin(), ip));

        ASSERT((OutsideIndices_.size() <= vertices.size()), "The outside indices size = " << OutsideIndices_.size() << " while vertices size = " << vertices.size());
    }
    else
    {
        OutsideIndices_.clear();
        OutsideIndices_.resize(vertices.size());
        std::iota(OutsideIndices_.begin(), OutsideIndices_.end(),0);
    }
}

void MeshCurvature::refine(Mesh& mesh)
{
    mesh_ = &mesh;

    // find the vertices to be updated
    findVertices();

    // Start refining 
    err_ = 10000000;
    iteration_=1;

    while (err_ >= tol_)
    {
        auto start = std::chrono::high_resolution_clock::now();

        Real maxerr = -100000;

        // first let's calculate the curvatures 
        curvatureCalc_ -> calculate(*mesh_);

        // then let's obtain the curvatures for each vertices 
        const auto& curvatures = curvatureCalc_ -> getAvgCurvaturePerVertex();

        // get vertics 
        auto& vertices = mesh_->accessvertices();

        // update the vertices 
        std::cout << "Iteration " << iteration_ << std::endl;
        Real avgE=0.0;
        #pragma omp parallel
        {
            Real e=maxerr;
            Real avgElocal = 0.0;
            #pragma omp for
            for (int j=0;j<OutsideIndices_.size();j++)
            {
                int index = OutsideIndices_[j];
                Real diffkappa = meanCurvature_ - curvatures[index];
 
                avgElocal += std::abs(diffkappa);

                if (std::abs(diffkappa) > e)
                {
                    e = std::abs(diffkappa);
                }

                for (int k=0;k<3;k++)
                {
                    vertices[index].position_[k] += StepSize_ * vertices[index].normals_[k] * diffkappa;
                }
            }

            #pragma omp critical
            if (e > maxerr)
            {
                maxerr = e;
            }

            #pragma omp critical
            {
                avgE += avgElocal;
            }
        }

        avgE = avgE / vertices.size();

        if (meanCurvature_ == 0)
        {
            err_ = 0;
        }
        else
        {
            err_ = maxerr / meanCurvature_;
        }

        std::cout << "Max error is " << err_ << std::endl;
        std::cout << "average error is " << avgE << "\n";

        // update the normals as well 
        mesh_->CalcVertexNormals();

        // update the iterations
        iteration_ ++;

        // break the calculation if iteration is already maxed out
        if (iteration_ > maxStep)
        {
            break;
        }

        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
        std::cout << "update = " << duration.count() << " us " << "\n";
    }
}