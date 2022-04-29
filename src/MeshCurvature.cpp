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
    // first let's see if user has provided with information about bounding boxes 
    Real2 xrange, yrange, zrange;
    Real INF = std::numeric_limits<Real>::infinity();

    xrange = {{-INF, INF}};
    yrange = {{-INF, INF}};
    zrange = {{-INF, INF}};

    bool xrangeRead = pack_.ReadArrayNumber("xrange", ParameterPack::KeyType::Optional, xrange);
    bool yrangeRead = pack_.ReadArrayNumber("yrange", ParameterPack::KeyType::Optional, yrange);
    bool zrangeRead = pack_.ReadArrayNumber("zrange", ParameterPack::KeyType::Optional, zrange);

    const auto& vertices = mesh_->getvertices();

    // If none of these are provided, then we keep the boundary vertices constant
    if (fixBoundary_)
    {
        mesh_->MapEdgeToFaces();
        mesh_->findBoundaryVertices();

        for (int i=0;i<vertices.size();i++)
        {
            if (! mesh_->isBoundary(i))
            {
                OutsideIndices_.push_back(i);
            }
        }
    }

    if (xrangeRead || yrangeRead || zrangeRead)
    {
        const auto& vertices = mesh_->getvertices();

        for (int i=0;i<vertices.size();i++)
        {
            auto v = vertices[i];
            bool inX = (v.position_[0] > xrange[0]) && (v.position_[0] < xrange[1]);
            bool inY = (v.position_[1] > yrange[0]) && (v.position_[1] < yrange[1]);
            bool inZ = (v.position_[2] > zrange[0]) && (v.position_[2] < zrange[1]);

            // if not all of them are inside the thing
            if (! (inX && inY && inZ))
            {
                OutsideIndices_.push_back(i);
            }
        }
    }

    // Find the unique indices 
    auto ip = std::unique(OutsideIndices_.begin(), OutsideIndices_.end());
    OutsideIndices_.resize(std::distance(OutsideIndices_.begin(), ip));

    ASSERT((OutsideIndices_.size() <= vertices.size()), "The outside indices size = " << OutsideIndices_.size() << " while vertices size = " << vertices.size());
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