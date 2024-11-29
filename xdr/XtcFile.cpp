#include "XtcFile.h"
namespace XdrFiles
{
    registry_<XtcFile> register_Xtc("xtc");
}

XtcFile::XtcFile(const XdrInput& input)
:XdrWrapper(input)
{}

XtcFile::XtcFile(std::string filename, std::string mode)
:XdrWrapper(filename, mode)
{

}


void XtcFile::readNumAtoms()
{
    ASSERT((isOpen()), "The file is not opened.");

    int success = read_xtc_natoms(const_cast<char*>(path_.c_str()), &natoms_);

    ASSERT((success == exdrOK), "The process to read xtc natoms is not sucessful.");
}

void XtcFile::writeFrame(const std::vector<Real3>& pos, int step, Real time, Matrix box)
{
    int size = pos.size();
    std::vector<rvec> tempRvec(size);

    for (int i=0;i<size;i++)
    {
        for (int j=0;j<3;j++)
        {
            tempRvec[i][j] = pos[i][j];
        }
    }

    matrix mat;
    for (int i=0;i<3;i++)
    {
        for (int j=0;j<3;j++)
        {
            mat[i][j] = box[i][j];
        }
    }

    int res = write_xtc(file_, size, step, time, mat, tempRvec.data(), 3);
    ASSERT((res == exdrOK), "Write failed.");
}


void XtcFile::readNframes()
{
    offsets_.clear();
    int est_nframes=0;
    int64_t* offsets=nullptr;

    read_xtc_n_frames(const_cast<char*>(path_.c_str()),&nframes_, &est_nframes, &offsets);
    offsets_.insert(offsets_.end(),offsets, offsets+nframes_);
    writeCheckPoint(path_, offsets_);
    std::cout << "nframes = " << nframes_ << std::endl;
}

void XtcFile::readFrame(int FrameNum)
{
    ASSERT((isOpen()), "The file is not opened.");
    ASSERT((FrameNum < nframes_), "Frame exceeded the maximum frames in xtc file");

    auto& positions_ = frame_.accessPositions();
    ASSERT((positions_.size() == natoms_), "The shape of positions is not of natoms size.");

    int step;
    Frame::Real time;
    matrix box_;
    Frame::Matrix matrix_box_;

    auto positions_ptr = (rvec*)positions_.data();

    // seek to the offset number 
    xdr_seek(file_,offsets_[FrameNum], 0);

    // read the xtc file
    int ret = read_xtc(file_, natoms_, &step, (float*)&time, box_, positions_ptr, (float*)&precision_);
    ASSERT((ret == exdrOK || ret == exdrENDOFFILE), "The reading operation in xtc file is not sucessful.");

    for(int i=0;i<3;i++)
    {
        for(int j=0;j<3;j++)
        {
            matrix_box_[i][j] = box_[i][j];
        }
    }
    frame_.setBox(matrix_box_);
    frame_.setTime(time);
    frame_.setStep(step);
}