#include "AverageField.h"
#include "LinAlgTools.h"

#include <string>
#include <vector>
#include <algorithm>
#include <iostream>

class AverageCurvature:public AverageField
{
    public:
        using index2 = std::array<int,2>;
        AverageCurvature(const DensityFieldInput& input);
        virtual ~AverageCurvature(){};

        virtual void finishCalculate() override;
        virtual void printFinalOutput() override;
        

        void printCurvature();
        void CalcTriangleAreaAndFacetNormals();
        void CalcVertexNormals();
        void CalcCurvature();

        // A function that obtains the neighbor lists as well as number of neighbors for a particular vertex
        void CalcNeighborsOfVertex();
    
    private:    
        std::vector<Real> curvature_;
        std::vector<Real> TriangleAreas_;
        std::vector<Real3> TriangleNormals_;
        std::vector<Real3> vertexNormals_;

        std::vector<int> Num_neighbors_;

        // same size as neighbor_indices_ but also specify which triangle the neighbor and itself belongs to 
        std::vector<std::vector<int>> vertex_triangle_indices_;
        std::vector<std::vector<int>> neighbor_indices_;
        std::vector<Real> neighborAreaTotlist_;

        std::string curvatureOutputName_;
        std::ofstream curvatureofs_;
};