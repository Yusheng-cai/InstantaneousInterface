#pragma once

#include "DensityField.h"
#include "marching_cubes.hpp"

#include <vector>
#include <string>

/*@brief Class that calculates the instantaneous Field as well as the instantaneous interface of the system  
*/
class InstantaneousField : public DensityField 
{
    public:
        InstantaneousField(const DensityFieldInput& input);

        virtual void calculate();
    private:
        int index_;
};