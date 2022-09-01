#include "Driver.h"

Driver::Driver(const ParameterPack& pack, const CommandLineArguments& cmd)
: pack_(const_cast<ParameterPack&>(pack)), cmd_(const_cast<CommandLineArguments&>(cmd))
{
    bool read = cmd_.readString("abspath", CommandLineArguments::Keys::Optional, abs_path_);
    if ( ! read)
    {
        abs_path_ = FileSystem::getCurrentPath();
    }

    // required xdr file input
    initializeXdrFile();

    // find gro file input
    initializeGroFile();

    // find all the atomgroup packs
    initializeAtomGroups();

    // Find the bounding box
    initializeBoundingBox();

    // initialize curvature 
    initializeCurvature();

    // find the density field
    initializeDensityField();

    // find the driver pack which is optional
    initializeDriver();
}

void Driver::initializeDriver()
{
    auto driverPack = pack_.findParamPack("driver", ParameterPack::KeyType::Optional);

    if (driverPack != nullptr)
    {
        driverPack -> Readbool("RandomSample", ParameterPack::KeyType::Optional, RandomSample_);
    }

    // if not bootstrap, then we read different stuff from driver pack 
    if (! RandomSample_)
    {
        if (driverPack != nullptr)
        {
            driverPack -> ReadNumber("startingframe", ParameterPack::KeyType::Optional, starting_frame_);
            driverPack -> ReadNumber("skip", ParameterPack::KeyType::Optional, skip_);
            ASSERT((skip_ > 0), "Skip must be an non-negative number, but provided skip = " << skip_);
            ASSERT((starting_frame_ >= 1), "The starting frame must be a number larger or equal to 1, but provided \
            startng frame = " << starting_frame_);

            // We do the 0 based counting inside the code
            starting_frame_ -= 1;
        }

        // arange the simulation frames from starting frame to total frame
        SimulationFrames_ = Algorithm::arange(starting_frame_, Totalframes_-1, skip_);
        simstate_.setTotalFramesToBeCalculated(SimulationFrames_.size());
    }
    else
    {
        driverPack -> ReadNumber("NumberRandomSample", ParameterPack::KeyType::Required, numRandomSample_);
        Algorithm::Permutation(Totalframes_, numRandomSample_, SimulationFrames_);
        simstate_.setTotalFramesToBeCalculated(numRandomSample_);
    }
}

void Driver::initializeDensityField()
{
    auto DensityPack = pack_.findParamPack("densityfield", ParameterPack::KeyType::Required);

    std::string densityType;
    DensityPack->ReadString("type", ParameterPack::KeyType::Required, densityType);

    DensityFieldInput input = {simstate_,const_cast<ParameterPack&>(*DensityPack), reg_};
    densityfield_ = DensityPtr(DensityFieldRegistry::Factory::instance().create(densityType, input));
}

void Driver::initializeBoundingBox()
{
    auto bbPack = pack_.findParamPacks("boundingbox", ParameterPack::KeyType::Required);

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

void Driver::initializeMeshRefinement()
{
    auto refinePack = pack_.findParamPacks("refine", ParameterPack::KeyType::Optional);

    for (int i=0;i<refinePack.size();i++)
    {
        std::string type;
        std::string name;

        auto& p = refinePack[i];

        p -> ReadString("type", ParameterPack::KeyType::Required, type);

        MeshRefineStrategyInput input = {const_cast<ParameterPack&>(*p)};
        meshRefineptr ptr = meshRefineptr(MeshRefineStrategyFactory::Factory::instance().create(type, input));

        reg_.registerMeshRefinement(ptr -> getName(), *ptr);
    }

}

void Driver::initializeCurvature()
{
    auto Cpack = pack_.findParamPacks("curvature", ParameterPack::KeyType::Optional);

    for (int i=0;i<Cpack.size();i++)
    {
        std::string curvatureType;
        std::string name;
        auto& p = Cpack[i];

        p -> ReadString("type", ParameterPack::KeyType::Required, curvatureType);

        CurvatureInput input = { const_cast<ParameterPack&>(*(p))};
        curveptr ptr = curveptr(CurvatureRegistry::Factory::instance().create(curvatureType, input));

        curvatures_.push_back(std::move(ptr)); 

        reg_.registerCurvature(curvatures_[i] -> getName(), *curvatures_[i]);
    }
}

void Driver::readFrameXdr(int FrameNum)
{
    xdrfile_ -> readFrame(FrameNum);

    // set the simulation box
    auto& simbox = xdrfile_->getSimulationBox(); 
    simstate_.setSimulationBox(simbox);
    simstate_.setStep(xdrfile_->getStep());
    simstate_.setTime(xdrfile_->getTime());
}

void Driver::update(int FrameNum)
{
    // read the xdr frame 
    readFrameXdr(FrameNum);

    // obtain the positions of the atoms from the xdr file (xtc, trr etc.)
    const auto& positions_ = xdrfile_ -> getPositions();

    // update the atomgroups
    for (int i=0;i<AtomGroupNames_.size();i++)
    {
        auto  agName = AtomGroupNames_[i];
        auto& ag     = simstate_.getAtomGroup(agName);

        ag.update(positions_, FrameNum); 
    }
}

void Driver::initializeXdrFile()
{
    auto xdrPack = pack_.findParamPack("xdrfile",ParameterPack::KeyType::Required);

    std::string xdrPath;
    xdrPack -> ReadString("path", ParameterPack::KeyType::Required, xdrPath);

    int found = xdrPath.find_last_of(".");
    ASSERT((found != std::string::npos), "The extension of the file is not specified.");

    std::string type = xdrPath.substr(found+1);

    XdrInput input = {const_cast<ParameterPack&>(*xdrPack), abs_path_};

    xdrfile_ = XdrPtr(XdrFiles::factory::instance().create(type, input));

    xdrfile_ -> open();

    Totalframes_ = xdrfile_ -> getNframes();
    simstate_.setTotalFrames(Totalframes_);
}

void Driver::initializeGroFile()
{
    auto groPack = pack_.findParamPack("grofile", ParameterPack::KeyType::Optional);

    if (groPack != nullptr)
    {
        groPack->ReadString("path", ParameterPack::KeyType::Required, groPath_);

        grofile_.Open(groPath_);
    }
}

void Driver::initializeAtomGroups()
{
    auto agPack = pack_.findParamPacks("atomgroup", ParameterPack::KeyType::Required);

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

void Driver::printOutputfileIfOnStep()
{
    densityfield_->printOutputIfOnStep();
}

void Driver::printFinalOutput()
{
    densityfield_->printFinalOutput();
}

void Driver::run()
{
    for (int i=0;i<SimulationFrames_.size();i++)
    {
        int ind = SimulationFrames_[i];

        std::cout << "Frame = " << ind << std::endl;

        // perform an update on the atom positions etc.
        update(ind);

        auto start = std::chrono::high_resolution_clock::now();

        // perform the calculations
        calculate();

        printOutputfileIfOnStep();
        auto end = std::chrono::high_resolution_clock::now();
        auto diff = std::chrono::duration_cast<std::chrono::milliseconds>(end-start);
        std::cout << "Time it took for calculate is " << diff.count() << " milliseconds." << std::endl;
    }

    // finish the calculation step
    finishCalculate();

    // print the final output
    printFinalOutput();
}