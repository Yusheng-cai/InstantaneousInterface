#pragma once
#include "CommonTypes.h"
#include "Assert.h"

#include <functional>
#include <vector>
#include <map>
#include <array>
#include <iostream>
#include <fstream>

class Output
{
    public:
        using outputFunc = std::function<void(std::string name)>;
        using perIteroutputFunc = std::function<void(std::ofstream& ofs)>;
        using perIterIndoutputFunc = std::function<void(std::string name)>;

        Output() = default;

        void registerOutputFunc(std::string name, outputFunc func);
        outputFunc& getOutputFuncByName(std::string name);

        void registerPerIterOutputFunc(std::string name, perIteroutputFunc func);
        perIteroutputFunc& getPerIterOutputFuncByName(std::string name);

        void registerPerIterIndOutputFunc(std::string name, perIterIndoutputFunc func);
        perIterIndoutputFunc& getPerIterIndOutputFuncByName(std::string name);

        const std::vector<std::string>& getOutputNames() const { return OutputNames_;}
        const std::vector<std::string>& getPerIterOutputNames() const {return PerIterOutputNames_;}
        const std::vector<std::string>& getPerIterIndOutputNames() const {return PerIterIndOutputNames_;}

    private:
        std::map<std::string, outputFunc> MapNameToOutputFunc_;
        std::map<std::string, perIteroutputFunc> MapNameToPerIterOutputFunc_;
        std::map<std::string, perIterIndoutputFunc> MapNameToPerIterIndOutputFunc_;
        std::vector<std::string> OutputNames_;
        std::vector<std::string> PerIterOutputNames_;
        std::vector<std::string> PerIterIndOutputNames_;
};