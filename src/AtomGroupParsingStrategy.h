#pragma once
#include "xdr/GroFile.h"
#include "tools/GenericFactory.h"

#include <vector>
#include <string>
#include <sstream>
#include <map>
#include <algorithm>
#include <iostream>

struct AtomGroupParsingInput
{
    GroFile& grofile_;
    std::vector<std::string> selection_str;
};

class AtomGroupParsingStrategy
{
    public:
        AtomGroupParsingStrategy(AtomGroupParsingInput& input):grofile_(input.grofile_), selection_str_(input.selection_str){};
        virtual ~AtomGroupParsingStrategy(){};

        virtual void Parse(std::vector<int>& indices) = 0;
        virtual void update(std::vector<int>& indices, int FrameNum){};

        void SortAndCheckNoDuplicate(std::vector<int>& indices);

        virtual int getMaxNumAtoms() {ASSERT((parsed), "Trying to get the max length before parsing the file."); return MaxNumAtoms_;}

    protected:
        GroFile& grofile_;
        std::vector<std::string>& selection_str_;
        std::vector<std::string> index_str_;
        int MaxNumAtoms_;
        bool parsed=false;
};

class AtomIndexParsing:public AtomGroupParsingStrategy
{
    public:
        AtomIndexParsing(AtomGroupParsingInput& input):AtomGroupParsingStrategy(input){};
        virtual ~AtomIndexParsing(){};

        virtual void Parse(std::vector<int>& indices);
};

class ResidueNumberParsing:public AtomGroupParsingStrategy
{
    public:
        ResidueNumberParsing(AtomGroupParsingInput& input):AtomGroupParsingStrategy(input){};
        virtual ~ResidueNumberParsing(){};

        virtual void Parse(std::vector<int>& indices);
};

class AtomTypeParsing:public AtomGroupParsingStrategy
{
    public:
        AtomTypeParsing(AtomGroupParsingInput& input):AtomGroupParsingStrategy(input){};
        virtual ~AtomTypeParsing(){};

        virtual void Parse(std::vector<int>& indices);
};

class ResidueNameParsing:public AtomGroupParsingStrategy
{
    public:
        ResidueNameParsing(AtomGroupParsingInput& input):AtomGroupParsingStrategy(input){};
        virtual ~ResidueNameParsing(){};

        virtual void Parse(std::vector<int>& indices);
};

class AtomIndexFile : public AtomGroupParsingStrategy
{
    public:
        AtomIndexFile(AtomGroupParsingInput& input) : AtomGroupParsingStrategy(input){};
        virtual ~AtomIndexFile(){};

        virtual void Parse(std::vector<int>& indices);
        virtual void update(){};

    private:    
        std::string fileName_;
};

class IndexFileParsing: public AtomGroupParsingStrategy
{
    public:
        IndexFileParsing(AtomGroupParsingInput& input);
        virtual ~IndexFileParsing(){};

        virtual void update(std::vector<int>& indices, int FrameNum);
        virtual void Parse(std::vector<int>& indices);

        const std::vector<std::vector<int>>& getIndexFileIndices() const { return Fileindices_;} 
        std::vector<std::vector<int>>& accessIndexFileIndices() { return Fileindices_;}

        const std::vector<int>& getFrame() const {return Frames_;}
        std::vector<int>& getFrame() {return Frames_;}

        int getNumFrames() const {return Frames_.size();}

        bool isOpen();

    private:
        std::string fileName_;

        std::ifstream ifs_;
        std::stringstream ss_;

        std::vector<std::vector<int>> Fileindices_;

        std::string comment_symbol = "#";

        int frame_count = 0;

        int totalFrames_;

        // correct for whether or not the indexing is 1 or 0 based --> assumes 1-based 
        int reduce_ = 1;

        // number to be skipped
        int skip_ = 0;

        std::vector<int> Frames_;
};

namespace AtomGroupParsingRegistry
{
    using Base = AtomGroupParsingStrategy;
    using Key  = std::string;

    using Factory = GenericFactory<Base,Key,AtomGroupParsingInput&>;

    template<typename D>
    using registry_ = RegisterInFactory<Base, D, Key, AtomGroupParsingInput&>;
}