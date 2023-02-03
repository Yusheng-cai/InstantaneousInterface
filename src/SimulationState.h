#pragma once
#include "tools/CommonTypes.h"
#include "SimulationBox.h"
#include "AtomGroup.h"
#include "BoundingBox.h"

#include <string>
#include <map>
#include <memory>

// This is a class that keeps track of the simulation progression
class SimulationState
{
    public:
        using Real = CommonTypes::Real;
        using Real3= CommonTypes::Real3;
        using Matrix=CommonTypes::Matrix;

        SimulationState(){};
        ~SimulationState(){};

        void setSimulationBox(Matrix boxMat){box_.setBoxMatrix(boxMat);}
        void setTime(Real time){time_ = time;}
        void setStep(int step){step_ = step;}
        void setFrameNum(int frame){frame_=frame;}
        void setTotalFrames(int totalFrames){TotalFrames_ = totalFrames;}
        void setTotalFramesToBeCalculated(int totalframestobecalculated){TotalFramesToBeCalculated_ = totalframestobecalculated;}

        // getters
        SimulationBox& getSimulationBox(){return box_;}
        const SimulationBox& getSimulationBox() const{return box_;}
        Real getTime() const{return time_;}
        int getStep() const {return step_;}
        int getFrame() const {return frame_;}
        int getTotalFrames() const {return TotalFrames_;}
        int getTotalFramesToBeCalculated() const {return TotalFramesToBeCalculated_;}

        // registers AtomGroups
        void registerAtomGroup(std::string name, AtomGroup& ag);

        // registers bounding box (for instantaneous interface)
        void registerBoundingBox(std::string name, BoundingBox& box);

        // get AtomGroup reference by name
        const AtomGroup& getAtomGroup(std::string name) const;
        AtomGroup& getAtomGroup(std::string name);

        // get bounding box by name
        const BoundingBox& getBoundingBox(std::string name) const;
        BoundingBox& getBoundingBox(std::string name);

        const std::map<std::string, AtomGroup>& getAtomGroupRegistry() const{ return MapName2AtomGroup_;}
        const std::map<std::string, BoundingBox>& getBoundingBoxRegistry() const {return MapName2BoundingBox_;}

    private:
        Real time_;
        int step_;
        int TotalFrames_;
        int frame_;
        int TotalFramesToBeCalculated_;

        SimulationBox box_;
        std::map<std::string,AtomGroup> MapName2AtomGroup_;
        std::map<std::string, BoundingBox> MapName2BoundingBox_;
};