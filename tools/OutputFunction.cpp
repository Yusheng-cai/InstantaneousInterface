#include "OutputFunction.h"

void Output::registerOutputFunc(std::string name, outputFunc func)
{
    auto it = MapNameToOutputFunc_.find(name);

    ASSERT((it == MapNameToOutputFunc_.end()), "The output with name " << name << " is already registered.");

    OutputNames_.push_back(name);

    MapNameToOutputFunc_.insert(std::make_pair(name, func));
}

Output::outputFunc& Output::getOutputFuncByName(std::string name)
{
    auto it = MapNameToOutputFunc_.find(name);

    ASSERT((it != MapNameToOutputFunc_.end()), "The output with name " << name << " is not registered.");

    return it -> second;
}

void Output::registerPerIterOutputFunc(std::string name, perIteroutputFunc func)
{
    auto it = MapNameToPerIterOutputFunc_.find(name);

    ASSERT((it == MapNameToPerIterOutputFunc_.end()), "The per iter output with name " << name << " is already registered.");

    PerIterOutputNames_.push_back(name);

    MapNameToPerIterOutputFunc_.insert(std::make_pair(name, func));
}

Output::perIteroutputFunc& Output::getPerIterOutputFuncByName(std::string name)
{
    auto it = MapNameToPerIterOutputFunc_.find(name);

    ASSERT((it != MapNameToPerIterOutputFunc_.end()), "The output with name " << name << " is not registered.");

    return it -> second;
}
void Output::registerPerIterIndOutputFunc(std::string name, perIterIndoutputFunc func)
{
    auto it = MapNameToPerIterIndOutputFunc_.find(name);

    ASSERT((it == MapNameToPerIterIndOutputFunc_.end()), "The per iter output with name " << name << " is already registered.");

    PerIterIndOutputNames_.push_back(name);

    MapNameToPerIterIndOutputFunc_.insert(std::make_pair(name, func));
}

Output::perIterIndoutputFunc& Output::getPerIterIndOutputFuncByName(std::string name)
{
    auto it = MapNameToPerIterIndOutputFunc_.find(name);

    ASSERT((it != MapNameToPerIterIndOutputFunc_.end()), "The output with name " << name << " is not registered.");

    return it -> second;
}