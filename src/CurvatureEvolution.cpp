#include "CurvatureEvolution.h"

namespace MeshRefineStrategyFactory
{
    registry_<CurvatureEvolution> registercurve("CurvatureEvolution");
}

CurvatureEvolution::CurvatureEvolution(MeshRefineStrategyInput& input)
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
}


void CurvatureEvolution::findVertices()
{
    const auto& vertices = mesh_->getvertices();

    MeshTools::MapEdgeToFace(*mesh_, MapEdgeToFace_, MapVertexToEdge_);
    MeshTools::CalculateBoundaryVertices(*mesh_, MapEdgeToFace_, boundaryIndicator_);

    for (int i=0;i<vertices.size();i++)
    {
        if (! MeshTools::IsBoundary(i, boundaryIndicator_))
        {
            VertexIndices_.push_back(i);
        }
    }
}

void CurvatureEvolution::refine(Mesh& mesh)
{
    mesh_ = &mesh;

    // find the vertices to be updated
    findVertices();

    // Start refining 
    err_ = 10000000;
    iteration_=1;

    while (err_ >= tol_)
    {
        // set max error to be some random negative number for now 
        Real maxerr = -1;

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
            for (int j=0;j<VertexIndices_.size();j++)
            {
                int index = VertexIndices_[j];

                // find the difference in curvature 
                Real diffkappa = meanCurvature_ - curvatures[index];

                // add to local average error
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
            std::cout << "Evolution finished premature at iteration " << iteration_ << "\n";
            break;
        }
    }
}