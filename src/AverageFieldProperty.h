#pragma once
#include "DensityField.h"
#include "parallel/OpenMP_buffer.h"
#include "Curvature.h"

#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <memory>

class AverageFieldProperty:public DensityField
{
    public:
        using CurvaturePtr = std::unique_ptr<Curvature>;

        AverageFieldProperty(const DensityFieldInput& input);
        virtual ~AverageFieldProperty(){};

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
        std::string property_filename_;
        std::vector<std::vector<Real>> properties_;

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