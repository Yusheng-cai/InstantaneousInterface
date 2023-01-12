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
        void average();

        // print the entire 3d field (xyz)
        void printField(std::string name);

        // print the 2d field (xy, xz or yz)
        void print2dField(std::string name);
    private:

        std::map<std::string, int> Map2dNameToDimension = \
        {
            {"xy", 2}, 
            {"yz", 0}, 
            {"xz", 1}, 
            {"yx", 2}, 
            {"zy", 0}, 
            {"zx", 1}     
        };
};