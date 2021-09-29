#include "OutputFunction.h"

void Output::registerOutputFunc(std::string name, outputFunc func)
{
    auto it = MapNameToOutputFunc_.find(name);

    ASSERT((it == MapNameToOutputFunc_.end()), "The output with name " << name << " is already registered.");

    MapNameToOutputFunc_.insert(std::make_pair(name, func));
}

Output::outputFunc& Output::getOutputFuncByName(std::string name)
{
    auto it = MapNameToOutputFunc_.find(name);

    ASSERT((it != MapNameToOutputFunc_.end()), "The output with name " << name << " is not registered.");

    return it -> second;
}