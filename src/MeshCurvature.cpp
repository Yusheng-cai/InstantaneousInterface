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

    // read in xtc info
    xtcWrite_ = pack_.ReadString("xtcfile", ParameterPack::KeyType::Optional, xtcName);

    if (xtcWrite_)
    {
        pack_.ReadNumber("xtcskip", ParameterPack::KeyType::Optional, skip_);
        xtcFile_ = xtcptr(new XtcFile(xtcName, "write"));
        xtcFile_->open();
    }
}

void MeshCurvature::refine()
{
    // First find the map from edges to faces 
    mesh_.MapEdgeToFaces();

    // From the above map, find boundary vertices 
    mesh_.findBoundaryVertices();

    // Start refining 
    Real err = 10000000;
    int iterations=0;

    while (err >= tol_)
    {
        Real maxerr = -100000;

        // first let's calculate the curvatures 
        curvatureCalc_ -> calculate(mesh_);

        // then let's obtain the curvatures for each vertices 
        const auto& curvatures = curvatureCalc_ -> getAvgCurvaturePerVertex();

        // get vertics 
        auto& vertices = mesh_.accessvertices();

        // update the vertices 
        std::cout << "Iteration " << iterations << std::endl;
        #pragma omp parallel
        {
            Real e=maxerr;
            #pragma omp for
            for (int j=0;j<vertices.size();j++)
            {
                if ( ! mesh_.isBoundary(j))
                {
                    Real diffkappa = meanCurvature_ - curvatures[j];

                    if (std::abs(diffkappa/curvatures[j]) > e)
                    {
                        e = std::abs(diffkappa/meanCurvature_);
                    }

                    for (int k=0;k<3;k++)
                    {
                        vertices[j].position_[k] += StepSize_ * vertices[j].normals_[k] * diffkappa;
                    }
                }
            }

            #pragma omp critical
            if (e > maxerr)
            {
                maxerr = e;
            }
        }

        std::cout << "Max error is " << maxerr << std::endl;
        err = maxerr;

        // update the normals as well 
        mesh_.CalcVertexNormals();

        std::vector<Real3> vertPos(vertices.size());
        for (int j=0;j<vertices.size();j++)
        {
            vertPos[j] = vertices[j].position_;
        }

        iterations ++;

        if (xtcWrite_)
        {
            if (iterations % skip_ == 0)
            {
                xtcFile_->writeFrame(vertPos, xtcstep_, iterations,box);

                xtcstep_ ++;
            }
        }

        if (iterations > maxStep)
        {
            break;
        }
    }
}