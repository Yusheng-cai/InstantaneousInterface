#pragma once
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
        virtual void printOutputIfOnStep() override {};

        void printField(std::string name);
    private:
};