#include "Registry.h"

void Registry::registerCurvature(std::string name, const Curvature& curve)
{
    auto it = MapNameToCurvature_.find(name);
    ASSERT((it == MapNameToCurvature_.end()), "The name for curvature " << name << " is already registered.");

    MapNameToCurvature_.insert(std::make_pair(name, const_cast<Curvature*>(&curve)));
}

Curvature& Registry::getCurvature(std::string name)
{
    auto it = MapNameToCurvature_.find(name);

    ASSERT((it != MapNameToCurvature_.end()), "The name for curvature " << name << " is not registered.");

    return *(it->second);
}