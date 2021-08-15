#include "DensityField.h"
#include "parallel/OpenMP_buffer.h"
#include "Curvature.h"

#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <memory>

class AverageField:public DensityField
{
    public:
        using CurvaturePtr = std::unique_ptr<Curvature>;

        AverageField(const DensityFieldInput& input);
        virtual ~AverageField(){};

        virtual void calculate() override;
        virtual void update() override{};
        virtual void finishCalculate() override;
        virtual void printOutputIfOnStep() override; 
        virtual void printFinalOutput() override;

        void initializeCurvature(std::vector<const ParameterPack*>& pack);

        void printField();
        void printTriangleIndices();
        void printVertices();
        void printNormals();


    private:
        std::string fieldOutputFileName_;
        std::string triangleIndicesFileName_;
        std::string vertexFileName_;
        std::string vertexNormalFileName_;

        std::ofstream fieldofs_;
        std::ofstream triangleIndicesofs_;
        std::ofstream vertexofs_;
        std::ofstream vertexNormalofs_;

        std::vector<CurvaturePtr> curvatures_;
};