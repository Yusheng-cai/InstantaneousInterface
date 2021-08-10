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

        void precalculateDensity();

    private:
        std::vector<int> AtomIndicesInside_;
        OpenMP::OpenMP_buffer<std::vector<int>> AtomIndicesBuffer_;
        std::string outputFileName_;

        std::vector<Real> DensityPreCalculate_;
};