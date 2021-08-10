#include "Driver.h"

Driver::Driver(const ParameterPack& pack, const CommandLineArguments& cmd)
{
    bool read = cmd.readString("abspath", CommandLineArguments::Keys::Optional, abs_path_);
    if ( ! read)
    {
        abs_path_ = FileSystem::getCurrentPath();
    }

    // required xdr file input
    auto xdrPack = pack.findParamPack("xdrfile",ParameterPack::KeyType::Required);
    initializeXdrFile(xdrPack);

    // find gro file input
    auto groPack = pack.findParamPack("grofile", ParameterPack::KeyType::Optional);
    if(groPack != nullptr)
    {
        initializeGroFile(groPack);
    }

    // find all the atomgroup packs
    auto atomgroupPack = pack.findParamPacks("atomgroup", ParameterPack::KeyType::Required);
    initializeAtomGroups(atomgroupPack);

    // Find the bounding box
    auto boundingboxPack = pack.findParamPacks("boundingbox", ParameterPack::KeyType::Required);
    initializeBoundingBox(boundingboxPack);

    // find the density field
    auto DensityPack = pack.findParamPack("densityfield", ParameterPack::KeyType::Required);
    initializeDensityField(DensityPack);
}

void Driver::initializeDensityField(const ParameterPack* densityPack)
{
    std::string densityType;
    densityPack->ReadString("type", ParameterPack::KeyType::Required, densityType);

    DensityFieldInput input = {simstate_,const_cast<ParameterPack&>(*densityPack)};
    densityfield_ = DensityPtr(DensityFieldRegistry::Factory::instance().create(densityType, input));
}

void Driver::initializeBoundingBox(std::vector<const ParameterPack*>& bbPack)
{
    for (int i=0;i<bbPack.size();i++)
    {
        std::string bbName;
        auto pack = bbPack[i];

        pack -> ReadString("name", ParameterPack::KeyType::Required, bbName);
        BoundingBoxInput input = {const_cast<ParameterPack&>(*pack), simstate_};

        BoundingBox box(input);

        simstate_.registerBoundingBox(bbName, box);
    }
}

void Driver::update()
{
    xdrfile_ -> readNextFrame();

    // set the simulation box
    auto& simbox = xdrfile_->getSimulationBox(); 
    simstate_.setSimulationBox(simbox);

    auto& positions_ = xdrfile_ -> getPositions();

    for (int i=0;i<AtomGroupNames_.size();i++)
    {
        auto& agName = AtomGroupNames_[i];
        auto& ag     = simstate_.getAtomGroup(agName);

        ag.update(positions_); 
    }
}

void Driver::initializeXdrFile(const ParameterPack* xdrPack)
{
    std::string xdrPath;
    xdrPack -> ReadString("path", ParameterPack::KeyType::Required, xdrPath);

    int found = xdrPath.find_last_of(".");
    ASSERT((found != std::string::npos), "Somehow the . character is not found in the path, it is needed to specify\
    what kind of file is being passed in.");

    std::string type = xdrPath.substr(found+1);

    XdrInput input = {const_cast<ParameterPack&>(*xdrPack), abs_path_};

    xdrfile_ = XdrPtr(XdrFiles::factory::instance().create(type, input));

    xdrfile_ -> open();
}

void Driver::initializeGroFile(const ParameterPack* groPack)
{
    groPack->ReadString("path", ParameterPack::KeyType::Required, groPath_);

    grofile_.Open(groPath_);
}

void Driver::initializeAtomGroups(std::vector<const ParameterPack*>& agPack)
{
    for (int i=0;i<agPack.size();i++)
    {
        auto pack = agPack[i];
        std::string atomName_;
        pack -> ReadString("name", ParameterPack::KeyType::Required, atomName_);

        AtomGroupInput input = {const_cast<ParameterPack&>(*pack), grofile_};

        AtomGroup group(input);

        simstate_.registerAtomGroup(atomName_, group);

        AtomGroupNames_.push_back(atomName_);
    }
}

const AtomGroup& Driver::getAtomGroup(const std::string& name) const 
{
    return simstate_.getAtomGroup(name);
}

AtomGroup& Driver::getAtomGroup(const std::string& name) 
{
    return simstate_.getAtomGroup(name);
}

void Driver::calculate()
{
    densityfield_->calculate();
}

void Driver::finishCalculate()
{
    densityfield_->finishCalculate();
}

void Driver::printOutputfile()
{
    densityfield_->printOutputIfOnStep();
}