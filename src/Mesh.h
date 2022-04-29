#pragma once
#include "tools/Assert.h"
#include "tools/CommonTypes.h"
#include "LinAlgTools.h"
#include "tools/InputParser.h"
#include "MeshRefineStrategy.h"
#include "tools/OutputFunction.h"
#include "happly.h"

#include <vector>
#include <array>
#include <cmath>
#include <algorithm>
#include <sstream>
#include <memory>
#include <unordered_map>
#include <iomanip>


struct vertex
{
    using Real = CommonTypes::Real;
    using Real3= CommonTypes::Real3;
    using index2 = CommonTypes::index2;

    Real3 position_;
    Real3 normals_;
    // The gaussian curvature of this particular vertex
    Real Gcurvature_=0.0;
    // The average curvature of this particular vertex 
    Real Acurvature_=0.0;
    int index;
};

struct edge
{
    using Real = CommonTypes::Real;
    using Real3= CommonTypes::Real3;

    // two vertices make up an edge
    vertex vertex1_;
    vertex vertex2_;
};

struct triangle
{
    using  INT3 = CommonTypes::index3;

    INT3 triangleindices_;
    
    // each triangle has 3 edges
    std::array<edge,3> edges_;

    // each triangle has 3 vertices
    std::array<vertex,3> vertices_;
};

class Mesh
{
    public:
        using Real3 = CommonTypes::Real3;
        using Real  = CommonTypes::Real;
        using refinePtr = std::unique_ptr<MeshRefineStrategy>;
        using INT2  = CommonTypes::index2;
        using INT3  = CommonTypes::index3;

        Mesh(const ParameterPack* pack);
        Mesh() {registerFunc();};
        ~Mesh(){};

        void registerFunc();

                                        /*******************************************
                                         ********** Printing Function **************
                                         ******************************************/

        void printSTL(std::string name);
        void printPLY(std::string name);
        void printPLYlibr(std::string name);
        void printPLYlibrCurvature(std::string name);
        void printBoundaryVertices(std::string name);
        void printArea(std::string name);
        void printCuttedMesh(std::string name);
        void printNonPBCMesh(std::string name);
        void printTranslatedMesh(std::string name);
        void printNeighbors(std::string name);

                                        /*******************************************
                                         * ******** setting Function **************
                                         * ****************************************/
        void setBoxLength(Real3 box) { boxLength_ = box; isPeriodic_=true;}

                                        /********************************************
                                         * ******** accessing Function **************
                                         * *****************************************/
        std::vector<vertex>& accessvertices() {return vertices_;}
        std::vector<triangle>& accesstriangles() {return triangles_;}
        std::vector<Real>& accessTriangleArea() {return triangleArea_;}
        std::vector<Real3>& accessPerVertexDir1() {return PerVertexdir1_;}
        std::vector<Real3>& accessPerVertexDir2() {return PerVertexdir2_;}
        Output& accessOutput() {return outputs_;}

        // getters
        const std::vector<vertex>& getvertices() const {return vertices_;}
        const std::vector<triangle>& gettriangles() const {return triangles_;}
        const std::vector<Real>& getTriangleArea() const {return triangleArea_;}
        const std::vector<Real3>& getPerVertexDir1() const {return PerVertexdir1_;}
        const std::vector<Real3>& getPerVertexDir2() const {return PerVertexdir2_;}
        const std::vector<Real3>& getFaceNormals() const {return facetNormals_;}
        const std::vector<std::vector<int>>& getMapVertexToFace() const {return MapVertexIndicesToFaceIndices_;}
        const std::vector<Real3>& getCornerAreas() const {return cornerArea_;}
        int getNumVertices() const {return vertices_.size();}
        int getNumTriangles() const {return triangles_.size();}
        Real3 getBoxLength() {return boxLength_;}
        const std::map<INT2, std::vector<int>> getMapEdgeToFace() const {return MapEdgeToFace_;}

        // Scale all the vertices by a single number 
        void scaleVertices(Real num);

        // function that updates the normals of a mesh
        void updateNormals();

        // Calculate the corner area according to trimesh code C++
        void CalculateCornerArea();

        // printoutput
        void print();

        // convert pbc mesh to non pbs mesh  --> usually for visualization purpose, so let's not print the normals 
        // we don't actually change the vertices in the mesh obj
        void ConvertToNonPBCMesh(std::vector<Real3>& vertices, std::vector<INT3>& faces);

        // calculate volume
        Real calculateVolume();

        // function that finds the area of all the triangules on the surface --> Now takes into account PBC
        void CalcTriangleAreaAndFacetNormals();

        // function that calculates the normals of each of the vertex --> weighted by area 
        void CalcVertexNormalsAreaWeighted();

        // function that calculates the normals of each of the vertex --> not weighted by anything and updates the normals in each vertices
        void CalcVertexNormals();

        // Find the vertex direction 
        void CalcPerVertexDir();

        // update the triangles/vertex/edges if they were to be changed
        void update();

        // Check whether or not the mesh is periodic 
        bool isPeriodic() {return isPeriodic_;}

        // find PBC distance between 2 vertices 
        void getVertexDistance(const vertex& v1, const vertex& v2, Real3& distVec, Real& dist);
        void getVertexDistance(const Real3& v1, const Real3& v2, Real3& distVec, Real& dist);
        void CalculateShift(const Real3& v1, const Real3& v2, Real3& shiftVec);

        // move a vertex into pbc box
        void MoveVertexIntoBox(const Real3& OldVertPos, Real3& NewVertexPos);

        // find shifted vertex position
        Real3 getShiftedVertexPosition(const vertex& v1, const vertex& v2);

        // find the shit that takes a position into the box 
        Real3 getShiftIntoBox(const Real3& v1);
    
    private:
        std::vector<vertex> vertices_;
        std::vector<triangle> triangles_;

        std::vector<std::vector<int>> VertexTriangleIndices_;

        std::vector<Real> triangleArea_;
        std::vector<Real3> cornerArea_;
        std::vector<Real3> facetNormals_;

        std::vector<Real3> vertexNormals_;

        std::vector<Real3> PerVertexdir1_;
        std::vector<Real3> PerVertexdir2_;

        // map the edge to face 
        std::map<INT2, std::vector<int>> MapEdgeToFace_;

        // a map from vertex indices to the face indices  
        std::vector<std::vector<int>> MapVertexIndicesToFaceIndices_;

        // output function
        Output outputs_;

        // outputs
        std::vector<std::string> outs_;

        // outputNames
        std::vector<std::string> outputNames_;

        // indices of the boundary vertices 
        std::vector<bool> MapBoundaryVertexIndicesToTrue_;

        // map vertex to edges 
        std::vector<std::vector<INT2>> MapVertexToEdges_;

        // norms of all the edges in the mesh
        std::vector<Real3> EdgeNorms_;

        // the factor to which all the vertices on the mesh will be multiplied by
        Real factor_=1.0;

        // parameter pack pointer
        const ParameterPack* pack_;

        // whether or not the mesh is periodic
        bool isPeriodic_=false;
        Real3 boxLength_;
};


namespace MeshTools
{
    using Real3 = CommonTypes::Real3;
    using Real  = CommonTypes::Real;
    using INT3  = CommonTypes::index3;
    using INT2  = CommonTypes::index2;

    bool readPLY(std::string& filename, Mesh& mesh_);
    bool readPLYlibr(std::string& filename, Mesh& mesh_);

    void writePLY(std::string filename, const std::vector<Real3>& Vertices, const std::vector<INT3>& faces, Real factor=1.0);
    void writePLY(std::string filename, const std::vector<Real3>& Vertices, const std::vector<INT3>& faces, const std::vector<Real3>& normals);

    // calculate PBC distance
    Real3 calculateShift(const Real3& vec1, const Real3& vec2, const Real3& boxLength);

    // find if an edge is periodic 
    bool isPeriodicEdge(const Real3& vec1, const Real3& vec2, Real3& newarr, const Real3& boxLength);

    // map vertices to faces
    void MapVerticesToFaces(Mesh& mesh, std::vector<std::vector<int>>& map);

    // calculate triangle areas and Face normals 
    void CalculateTriangleAreasAndFaceNormals(Mesh& mesh, std::vector<Real>& Areas, std::vector<Real3>& Normals);

    // find vertex neighbors 
    void CalculateVertexNeighbors(Mesh& mesh, std::vector<std::vector<int>>& neighborIndices);

    // map Edge to faces
    void MapEdgeToFace(Mesh& mesh, std::map<INT2,std::vector<int>>& magEdgeToFace, std::vector<std::vector<INT2>>& MapVertexToEdge);

    // find boundary vertices 
    void CalculateBoundaryVertices(Mesh& mesh, std::map<INT2, std::vector<int>>& mapEdgeToFace, std::vector<bool>& boundaryIndicator);

    // check if point is on boundary
    bool IsBoundary(int Index, const std::vector<bool>& boundaryIndicator);
};