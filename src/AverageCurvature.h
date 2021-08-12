#include "AverageField.h"
#include "LinAlgTools.h"

#include <string>
#include <vector>
#include <algorithm>
#include <iostream>

class AverageCurvature:public AverageField
{
    public:
        AverageCurvature(const DensityFieldInput& input);
        virtual ~AverageCurvature(){};

        virtual void finishCalculate() override;
        virtual void printFinalOutput() override;
        

        void printCurvature();
        void CalcTriangleAreaAndFacetNormals();
        void CalcVertexNormals();
        void CalcCurvature();
    
    private:    
        std::vector<Real> curvature_;
        std::vector<Real> TriangleAreas_;
        std::vector<Real3> TriangleNormals_;
        std::vector<Real3> vertexNormals_;

        std::vector<int> Num_neighbors_;

        std::string curvatureOutputName_;
        std::ofstream curvatureofs_;
};