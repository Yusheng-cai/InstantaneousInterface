#include "AtomGroup.h"
#include "SimulationState.h"
#include "DensityField.h"
#include "xdr/XdrWrapper.h"
#include "xdr/GroFile.h"
#include "tools/InputParser.h"
#include "tools/CommandLineArguments.h"
#include "tools/FileSystem.h"
#include "BoundingBox.h"

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

        Driver(const ParameterPack& pack, const CommandLineArguments& cmd);

        void initializeAtomGroups(std::vector<const ParameterPack*>& agPack);

        void initializeGroFile(const ParameterPack* groPack);

        void initializeXdrFile(const ParameterPack* xdrPack);

        void initializeBoundingBox(std::vector<const ParameterPack*>& bbPack);

        void initializeDensityField(const ParameterPack* densityPack);

        void update();

        void calculate();

        // get num frames
        int getNumFrames() const {return xdrfile_->getNframes();}

        const std::vector<Real3> getPositions() const {return xdrfile_->getPositions();}
        std::vector<Real3> accessPositions() {return xdrfile_->getPositions();}

        const AtomGroup& getAtomGroup(const std::string& name) const; 
        AtomGroup& getAtomGroup(const std::string& name);

    private:
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
};