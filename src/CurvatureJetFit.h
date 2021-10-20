#include "Curvature.h"
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Monge_via_jet_fitting.h>
#include "tools/CommonTypes.h"
#include "Graph.h"
#define CGAL_EIGEN3_ENABLED 1

#include <vector>

class CurvatureJetFit : public Curvature
{
    public:
        using DataKernel = CGAL::Simple_cartesian<double>;
        using Dpoint     = DataKernel::Point_3;
        using MongeViaJetFitting = CGAL::Monge_via_jet_fitting<DataKernel>;
        using Mongeform =  MongeViaJetFitting::Monge_form;
        using DVector    = DataKernel::Vector_3;

        CurvatureJetFit(CurvatureInput& input);
        virtual void calculate() override;

        virtual void printCurvature(std::string name) override;
        void printNeighbors(std::string name);
        void printCoefficientPerVertex(std::string name);

    private:
        MongeViaJetFitting jetfitter_;
        Mongeform mform_;
        int numneighbors_;
        int degree_=2;
        int MongeCoefficient_=2;

        std::vector<Real2> CurvaturePerVertex_;
        std::vector<std::vector<int>> NeighborIndicesNVertex_;
        std::vector<std::vector<Real>> coefficientsPerVertex_;
};
