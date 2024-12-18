#pragma once

#include "libxdr/xdrfile_trr.h"
#include "libxdr/xdrfile_xtc.h"
#include "libxdr/xdrfile.h"
#include "tools/Assert.h"
#include "Frame.h"
#include "tools/GenericFactory.h"
#include "tools/CommonTypes.h"
#include "tools/InputParser.h"
#include "GroFile.h"
#include "tools/FileSystem.h"
#include "tools/CommonOperations.h"

#include <string>
#include <vector>
#include <iostream>
#include <sstream>


struct XdrInput
{
    const ParameterPack& pack;
    // the absolute path
    std::string apath_;
};

class XdrWrapper
{
    public:
        using Real = CommonTypes::Real;
        using Real3= CommonTypes::Real3;
        using Matrix = CommonTypes::Matrix;

        enum Mode
        {
            Read, Write, Append
        };

        XdrWrapper(const XdrInput&);    
        XdrWrapper(std::string filename, std::string mode);
        void open();
        virtual ~XdrWrapper();

        // close the file
        void close();

        // Check if the xdr file is opened
        bool isOpen();

        int getNumAtoms() const {return natoms_;}

        // write check point file
        void writeCheckPoint(std::string path, std::vector<int64_t>& offsets);
        void readCheckPoint(std::string filename, std::vector<int64_t>& offsets);
        std::string CheckPointFileName(std::string filename);

        virtual void readFrame(int FrameNum) = 0;
        virtual void writeFrame(const std::vector<Real3>& pos, int step, Real time, Matrix box){};
        virtual void readNumAtoms()  = 0;
        virtual void readNframes(){};
        const Frame::VectorReal3& getPositions() const{return frame_.getPositions();}
        const Frame::VectorReal3& getVelocities() const{return frame_.getVelocities();}
        const Frame::VectorReal3& getForces() const{return frame_.getForces();}
        Real getTime() const {return frame_.getTime();}
        int getStep() const {return frame_.getStep();}
        int getNframes() const {return nframes_;}
        const Matrix& getSimulationBox() const {return frame_.getBoxMatrix();}
        
    protected:
        XDRFILE* file_=nullptr;
        int natoms_;
        std::string path_;
        std::string operation_mode_;

        // gro file information
        GroFile grofile_;
        std::string groname_="";

        Frame frame_;
        int nframes_=0;

        std::string apath_;

        // the offsets in terms of the number of bits
        std::vector<int64_t> offsets_;
};

namespace XdrFiles
{
    using Key  = std::string;
    using Base = XdrWrapper;

    using factory = GenericFactory<Base, Key, const XdrInput&>;

    template<class D>
    using registry_ = RegisterInFactory<Base, D, Key, const XdrInput&>;
};