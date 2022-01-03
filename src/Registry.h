#pragma once 

#include "tools/Assert.h"
#include "Curvature.h"

#include <string>
#include <map>
#include <vector>

class Registry
{
    public:
        Registry()  =default;

        void registerCurvature(std::string name, const Curvature& curve);
        Curvature& getCurvature(std::string name);

    private:
        std::map<std::string, Curvature*> MapNameToCurvature_;
};