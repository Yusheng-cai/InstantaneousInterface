#pragma once

#include "AtomGroup.h"
#include "SimulationState.h"
#include "DensityField.h"
#include "xdr/XdrWrapper.h"
#include "xdr/GroFile.h"
#include "tools/InputParser.h"
#include "tools/CommandLineArguments.h"
#include "tools/FileSystem.h"
#include "BoundingBox.h"
#include "Curvature.h"
#include "Registry.h"

#include <vector>
#include <memory>
#include <array>

class Driver
{
    public:
        using XdrPtr = std::unique_ptr<XdrWrapper>;
        using DensityPtr = std::unique_ptr<DensityField>;
        using BBPtr  = std::unique_ptr<BoundingBox>;
        using Real3  = CommonTypes::Real3;
        using curveptr = std::unique_ptr<Curvature>;
        using meshRefineptr = std::unique_ptr<MeshRefineStrategy>;

        Driver(const ParameterPack& pack, const CommandLineArguments& cmd);

        void initializeAtomGroups(std::vector<const ParameterPack*>& agPack);

        void initializeGroFile(const ParameterPack* groPack);

        void initializeXdrFile(const ParameterPack* xdrPack);

        void initializeBoundingBox(std::vector<const ParameterPack*>& bbPack);

        void initializeDensityField();

        void initializeDriver(const ParameterPack* driverPack);

        void initializeCurvature();

        // we initialize the parameters packs for mesh refinement
        void initializeMeshRefinement();

        // Read information from xdr
        void readFrameXdr(int FrameNum);

        void update();

        void calculate();

        void printOutputfileIfOnStep();

        void printFinalOutput();

        void finishCalculate();

        // check if the step is valid, the FrameNum is inputted in 0 based counting
        bool CheckValidStep(int FrameNum);

        // get num frames
        int getNumFrames() const {return xdrfile_->getNframes();}

        const std::vector<Real3> getPositions() const {return xdrfile_->getPositions();}
        std::vector<Real3> accessPositions() {return xdrfile_->getPositions();}

        const AtomGroup& getAtomGroup(const std::string& name) const; 
        AtomGroup& getAtomGroup(const std::string& name);

    private:
        ParameterPack pack_;

        SimulationState simstate_;
        GroFile grofile_;

        std::vector<std::string> AtomGroupNames_;
        std::vector<std::string> BBNames_;

        std::string groPath_;

        XdrPtr xdrfile_;

        // the absolute path of the program
        std::string abs_path_;

        DensityPtr densityfield_;

        BBPtr boundingbox_;

        // vector of curvature objects
        std::vector<curveptr> curvatures_;

        Registry reg_;

        int Totalframes_;

        // The starting frame_ + how many frames to skip for calculation
        // This is 1 based counting
        int starting_frame_ = 1;
        // the number of frames to be skipped between 2 calculation
        // if skip = 2, then if the first frame = 1, then the next frame to be calculate is 3
        int skip_ = 0;
};