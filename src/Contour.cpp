#include "Contour.h"

std::vector<Contour::Real3> Contour::ExtendContourLine(const std::vector<Real3>& ContourLine, Real extend){
    // find centroid
    Real3 centroid = {{0,0,0}};

    for (auto point : ContourLine){
        centroid = centroid + point;
    }

    centroid = centroid / ContourLine.size();

    // find the vector to centroid
    std::vector<Real3> vec_to_centroid(ContourLine.size());

    for (int i=0;i<ContourLine.size();i++){
        Real3 dist = ContourLine[i] - centroid;
        LinAlg3x3::normalize(dist);

        vec_to_centroid[i] = dist;
    }

    // extend the contour line
    std::vector<Real3> newContour(ContourLine.size());
    for (int i=0;i<ContourLine.size();i++){
        Real3 newP = ContourLine[i] + vec_to_centroid[i] * extend;
        newContour[i] = newP;
    }

    return newContour;
}