#pragma once
#include "CommonTypes.h"
#include "Assert.h"

#include <functional>
#include <vector>
#include <map>
#include <array>

class Output
{
    public:
        using outputFunc = std::function<void(std::string name)>;

        Output() = default;

        void registerOutputFunc(std::string name, outputFunc func);
        outputFunc& getOutputFuncByName(std::string name);

    private:
        std::map<std::string, outputFunc> MapNameToOutputFunc_;
};