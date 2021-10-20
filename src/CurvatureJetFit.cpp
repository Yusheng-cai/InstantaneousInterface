#include "CurvatureJetFit.h"

namespace CurvatureRegistry
{
    registry<CurvatureJetFit> registerJet("jetfit");
}

CurvatureJetFit::CurvatureJetFit(CurvatureInput& input)
:Curvature(input)
{
}

void CurvatureJetFit::calculate()
{
    auto& vertices = mesh_.accessvertices();
    std::vector<Dpoint> vec_;
    for (int i=0;i<vertices.size();i++)
    {
        auto& pos = vertices[i].position_;
        Dpoint point(pos[0], pos[1], pos[2]);
        vec_.push_back(point);
    }

    mform_ = jetfitter_(vec_.begin(), vec_.end(), 3, 2);

    for (int i=0;i<vertices.size();i++)
    {
        std::cout << "Curvature of " << i << " = " << mform_.principal_curvatures(i) << std::endl;
    }
    std::cout << "FINISHED." << std::endl;
}