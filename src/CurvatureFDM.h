#include "Curvature.h"

class CurvatureFDM: public Curvature
{
    public:
        CurvatureFDM(CurvatureInput& input);

        virtual void calculate();
        virtual ~CurvatureFDM(){};
        void calculateArithmeticMean();
        void calculateGeometricMean(); 

        virtual void printOutput();

    private:
        std::string MeanMethod_="geometric";
};