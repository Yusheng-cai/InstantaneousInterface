#include "Curvature.h"
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Monge_via_jet_fitting.h>
#include "tools/CommonTypes.h"
#define CGAL_EIGEN3_ENABLED 1

#include <vector>

class CurvatureJetFit : public Curvature
{
    public:
        using DataKernel = CGAL::Simple_cartesian<double>;
        using Dpoint     = DataKernel::Point_3;
        using MongeViaJetFitting = CGAL::Monge_via_jet_fitting<DataKernel>;
        using Mongeform =  MongeViaJetFitting::Monge_form;


        CurvatureJetFit(CurvatureInput& input);
        virtual void calculate() override;
    private:
        MongeViaJetFitting jetfitter_;
        Mongeform mform_;
};
