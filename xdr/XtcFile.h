#pragma once

#include "XdrWrapper.h"
#include "libxdr/xtc_seek.h"
#include "tools/FileSystem.h"

#include <chrono>
#include <memory>

class XtcFile:public XdrWrapper
{
    public:
        XtcFile(const XdrInput& input);
        XtcFile(std::string filename, std::string mode);

        virtual ~XtcFile(){};

        virtual void readFrame(int FrameNum) override;
        virtual void readNumAtoms() override;
        virtual void readNframes() override;
        virtual void writeFrame(const std::vector<Real3>& pos, int step, Real time, Matrix box);

    private:
        Frame::Real precision_;
};