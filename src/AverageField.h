#include "DensityField.h"
#include "parallel/OpenMP_buffer.h"

#include <string>
#include <vector>
#include <iostream>
#include <sstream>

class AverageField:public DensityField
{
    public:
        AverageField(const DensityFieldInput& input);

        virtual void calculate() override;
        virtual void update() {};
        virtual void finishCalculate() override;
        virtual void printOutputIfOnStep() override; 
        virtual void printFinalOutput() override;

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
};