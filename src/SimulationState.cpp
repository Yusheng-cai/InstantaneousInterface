#include "SimulationState.h"

void SimulationState::registerAtomGroup(std::string name, AtomGroup& ag)
{
    //MapName2AtomGroup_.insert(std::make_pair(name,std::move(ag)));
    MapName2AtomGroup_.emplace(name, std::move(ag));
}

void SimulationState::registerBoundingBox(std::string name, BoundingBox& box)
{
    MapName2BoundingBox_.emplace(name, std::move(box));
}

const AtomGroup& SimulationState::getAtomGroup(std::string name) const
{
    auto it = MapName2AtomGroup_.find(name);

    ASSERT((it != MapName2AtomGroup_.end()), "The name of AtomGroup: " << name << " is not registered in SimulationState.");

    return it -> second;
}

AtomGroup& SimulationState::getAtomGroup(std::string name)
{
    auto it = MapName2AtomGroup_.find(name);

    ASSERT((it != MapName2AtomGroup_.end()), "The name of AtomGroup: " << name << " is not registered in SimulationState.");

    return it -> second;
}

const BoundingBox& SimulationState::getBoundingBox(std::string name) const 
{
    auto it = MapName2BoundingBox_.find(name);

    ASSERT(( it != MapName2BoundingBox_.end()), "The name of Bounding Box: " << name << " is not registered in SimulationState.");

    return it -> second;
}

BoundingBox& SimulationState::getBoundingBox(std::string name) 
{
    auto it = MapName2BoundingBox_.find(name);

    ASSERT(( it != MapName2BoundingBox_.end()), "The name of bounding box: " << name << " is not registered in SimulationState.");

    return it -> second;
}