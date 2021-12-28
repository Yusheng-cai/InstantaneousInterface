#include "Curvature.h"

class CurvatureFDM: public Curvature
{
    public:
        CurvatureFDM(CurvatureInput& input);

        virtual void calculate(Mesh& mesh);
        virtual ~CurvatureFDM(){};
        void calculateArithmeticMean(Mesh& mesh);
        void calculateGeometricMean(Mesh& mesh); 

        virtual void printCurvature(std::string name) override;

    private:
        std::string MeanMethod_="geometric";
        std::vector<Real> curvature_;
};