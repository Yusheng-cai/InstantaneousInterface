#include "XdrWrapper.h"
XdrWrapper::XdrWrapper(const XdrInput& input)
: apath_(input.apath_)
{
    input.pack.ReadString("path", ParameterPack::KeyType::Required, path_);
    path_ = FileSystem::joinPath(apath_, path_);

    bool readmode = input.pack.ReadString("mode", ParameterPack::KeyType::Optional, operation_mode_);
 

    if ( ! readmode)
    {
        operation_mode_ = "read";
    }
}

XdrWrapper::XdrWrapper(std::string filename, std::string mode)
:operation_mode_(mode) , path_(filename)
{
}

void XdrWrapper::writeCheckPoint(std::string filename, std::vector<int>& offsets)
{
    std::string newOffsetName = CheckPointFileName(filename);

    std::ofstream ofs;
    ofs.open(newOffsetName);
    for (int i=0;i<offsets.size();i++)
    {
        ofs << offsets[i] << " ";
    }

    ofs.close();
}

std::string XdrWrapper::CheckPointFileName(std::string filename)
{
    std::string::size_type pos = filename.find(".");
    ASSERT((pos != std::string::npos), "The filename " << filename << " does not contain .");
    std::string xdrName = filename.substr(0, pos);
    std::string newOffsetName = "." + xdrName + "_offset.out";

    return newOffsetName;
}

void XdrWrapper::open()
{
    std::string mode_;

    ASSERT((operation_mode_ == "read" || operation_mode_ == "append" || operation_mode_ == "write"), "mode has to be one of 'read' & 'append' & 'write'");

    // default is reading mode
    if (operation_mode_ == "read")
    {
        mode_ = "r";
    }
    else if (operation_mode_ == "append")
    {
        mode_ = "a";
    }   
    else if (operation_mode_ == "write")
    {
        mode_ = "w";
    }

    file_ = xdrfile_open(path_.c_str(), mode_.c_str());
    ASSERT((isOpen()), "The file did not open correctly.");

    // Read Number of Atoms 
    if (mode_ == "r")
    {
        readNumAtoms();

        frame_.setNumAtoms(natoms_);

        readNframes();
    }
}


bool XdrWrapper::isOpen()
{
    if (file_ != nullptr)
    {
        return true;
    }
    else
    {
        return false;
    }
}

void XdrWrapper::close()
{
    ASSERT((isOpen()), "You are trying to close the file while it wasn't opened in the first place.");
    int ret = xdrfile_close(file_);

    ASSERT((ret == 0), "The xdr file was closed incorrectly.");
}

XdrWrapper::~XdrWrapper()
{
    close();
}