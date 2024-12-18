#include "TrrFile.h"

namespace XdrFiles
{
    static const registry_<TrrFile> register_trrfile("trr");
}

TrrFile::TrrFile(const XdrInput& input)
:XdrWrapper(input)
{};

void TrrFile::readNframes()
{
    offsets_.clear();

    // find the check point file name 
    std::string cptFile = CheckPointFileName(path_);

    // check if check point file exist
    if (FileSystem::FileExist(cptFile))
    {
        readCheckPoint(path_, offsets_);
        nframes_ = offsets_.size();
    }
    else
    {
        int est_nframes=0;
        int64_t* offsets=nullptr;

        read_trr_n_frames(const_cast<char*>(path_.c_str()),&nframes_, &est_nframes, &offsets);
        offsets_.insert(offsets_.end(),offsets, offsets+nframes_);
        writeCheckPoint(path_, offsets_);
    }
}

void TrrFile::readFrame(int FrameNum)
{
    ASSERT((isOpen()), "The file is not opened.");

    auto& position = frame_.accessPositions();
    ASSERT((position.size() == natoms_), "The position vector does not have the correct number of atoms.");
    auto& velocities = frame_.accessVelocities();
    ASSERT((velocities.size() == natoms_), "The velocity vector does not have the correct number of atoms.");
    auto& forces = frame_.accessForces();
    ASSERT((forces.size() == natoms_), "The force vector does not have the correct number of atoms.");

    rvec* position_ptr = (rvec*)position.data();
    rvec* velocities_ptr = (rvec*)velocities.data(); 
    rvec* forces_ptr = (rvec*)forces.data();

    matrix box;
    int step;
    Frame::Real time, lambda; 
    int has_prop;

    xdr_seek(file_, offsets_[FrameNum], 0);
    int ret = read_trr(file_, natoms_, &step, (float*)&time, (float*)&lambda, box, position_ptr, velocities_ptr, forces_ptr, &has_prop);

    frame_.setTime(time);
    frame_.setStep(step);
    Frame::Matrix box_;
    for (int i=0;i<3;i++)
    {
        for(int j=0;j<3;j++)
        {
            box_[i][j] = box[i][j];
        }
    }
    frame_.setBox(box_);
}

void TrrFile::readNumAtoms()
{
    ASSERT((isOpen()), "The file is not open.");

    int sucess = read_trr_natoms(const_cast<char*>(path_.c_str()),&natoms_);
    ASSERT((sucess == exdrOK), "Reading natoms failed.");
}