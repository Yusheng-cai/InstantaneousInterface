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

class MarchingCubesWrapper;


struct vertex
{
    using Real = CommonTypes::Real;
    using Real3= CommonTypes::Real3;

    Real3 position_;
    Real3 normals_;
    // The gaussian curvature of this particular vertex
    Real Gcurvature_=0.0;
    // The average curvature of this particular vertex 
    Real Acurvature_=0.0;
    int index;

    bool operator==(const vertex& v1) const
    {
        if (v1.index == index)
        {
            return true;
        }
        else
        {
            return false;
        }
    }
};

struct edge
{
    using Real = CommonTypes::Real;
    using Real3= CommonTypes::Real3;

    // two vertices make up an edge
    vertex vertex1_;
    vertex vertex2_;

    bool operator==(const edge& e1) const
    {
        // either vertex 1 == vertex 1 && vertex2 == vertex2
        if (e1.vertex1_ == vertex1_ && e1.vertex2_ == vertex2_)
        {
            return true;
        }
        // or vertex1 == vertex2 && vertex2 == vertex1
        if (e1.vertex1_ == vertex2_ && e1.vertex2_ == vertex1_)
        {
            return true;
        }
        else
        {
            return false;
        }
    }
};

namespace std {

  template <>
  struct hash<vertex>
  {
      std::size_t operator()(const vertex& v) const 
      {
          std::size_t indexHash = std::hash<int>()(v.index);

          return indexHash;
      }
  };

  template <>
  struct hash<edge>
  {
    std::size_t operator()(const edge& e) const
    {
      using std::size_t;
      using std::hash;
      using std::string;

      // Compute individual hash values for first,
      // second and third and combine them using XOR
      // and bit shifting:

      std::size_t v1_hash = std::hash<vertex>()(e.vertex1_);
      std::size_t v2_hash = std::hash<vertex>()(e.vertex2_);

      return v1_hash ^ v2_hash;
    };
  };
}

struct triangle
{
    using index3 = CommonTypes::index3;

    index3 triangleindices_;
    
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

        Mesh(const ParameterPack* pack);
        Mesh() {};
        ~Mesh(){};

        void initializeRefineStrat();

        void printSTL(std::string name);
        void printPLY(std::string name);
        void printPLYAng(std::string name);
        void printPLYnm(std::string name);
        void printPLYlibr(std::string name);
        void printPLYlibrCurvature(std::string name);

        std::vector<vertex>& accessvertices() {return vertices_;}
        std::vector<triangle>& accesstriangles() {return triangles_;}
        std::vector<int>& accessNumNeighbors() {return NumNeighbors_;}
        std::vector<std::vector<int>>& accessNeighborIndices() {return neighborIndices_;}
        std::vector<Real>& accessTriangleArea() {return triangleArea_;}
        std::vector<Real3>& accessPerVertexDir1() {return PerVertexdir1_;}
        std::vector<Real3>& accessPerVertexDir2() {return PerVertexdir2_;}

        const std::vector<vertex>& getvertices() const {return vertices_;}
        const std::vector<triangle>& gettriangles() const {return triangles_;}
        const std::vector<int>& getNumNeighbors() const {return NumNeighbors_;}
        const std::vector<std::vector<int>>& getNeighborIndices() const {return neighborIndices_;}
        const std::vector<Real>& getTriangleArea() const {return triangleArea_;}
        const std::vector<Real3>& getPerVertexDir1() const {return PerVertexdir1_;}
        const std::vector<Real3>& getPerVertexDir2() const {return PerVertexdir2_;}
        const std::unordered_map<vertex, std::vector<edge>>& getMapBVertexToBEdges() const {return MapBoundaryVertexToBoundaryEdges_;}
        const std::vector<Real3>& getFaceNormals() const {return facetNormals_;}
        const std::vector<std::vector<int>>& getMapVertexToFace() const {return MapVertexIndicesToFaceIndices_;}

        const std::vector<Real3>& getCornerAreas() const {return cornerArea_;}

        int getNumVertices() const {return vertices_.size();}
        int getNumTriangles() const {return triangles_.size();}
        bool isBoundary(int i);
        const std::unordered_map<edge, std::vector<int>> getMapEdgeToFace() const {return MapEdgeToFace_;}

        // Find neighbors indices for a vertex
        void findVertexNeighbors();

        // Scale all the vertices by a single number 
        void scaleVertices(Real num);

        // find all the boundary vertices
        void findBoundaryVertices();

        // function called when trying to refine a mesh
        void refine();

        // function that updates the normals of a mesh
        void updateNormals();

        void test();

        // printoutput
        void print();

        // calculate volume
        Real calculateVolume();

        // function that finds the triangle indices that a vertex belongs to 
        void findTriangleIndices();

        // function that finds the area of all the triangules on the surface
        void CalcTriangleAreaAndFacetNormals();

        // function that calculates the normals of each of the vertex
        void CalcVertexNormals();

        // Find the 
        void CalcPerVertexDir();

        // calculate the faces each edge corresponds to
        void MapEdgeToFaces();

        // map vertex to faces
        void MapVertexToFaces();

        // update the triangles/vertex/edges if they were to be changed
        void update();

        // get the edges corresponds to the particular vertex
        std::vector<edge>& getEdgeForVertex(int i);
        std::vector<int>& getFaceIndicesForEdge(const edge& e);
    
    private:
        std::vector<vertex> vertices_;
        std::vector<triangle> triangles_;

        // finds the neighbors of a vertex as well as count its number of neighbors
        std::vector<std::vector<int>> neighborIndices_;
        std::vector<int> NumNeighbors_;

        std::vector<std::vector<int>> VertexTriangleIndices_;

        std::vector<Real> triangleArea_;
        std::vector<Real3> cornerArea_;
        std::vector<Real3> facetNormals_;

        std::vector<Real3> vertexNormals_;

        std::vector<Real3> PerVertexdir1_;
        std::vector<Real3> PerVertexdir2_;

        std::unordered_map<edge, std::vector<int>> MapEdgeToFace_;

        std::string refineStrategy_;

        // ptr to the refine object
        refinePtr MeshRefine_;

        // a map from boundary vertices to their respective edges  
        std::unordered_map<vertex, std::vector<edge>> MapBoundaryVertexToBoundaryEdges_;

        // a map from vertex indices to the face indices  
        std::vector<std::vector<int>> MapVertexIndicesToFaceIndices_;

        // output function
        Output outputs_;

        // outputs
        std::vector<std::string> outs_;

        // outputNames
        std::vector<std::string> outputNames_;

        // indices of the boundary vertices 
        std::unordered_map<int,bool> MapBoundaryVertexIndicesToTrue_;

        // map vertex to edges 
        std::unordered_map<int, std::vector<edge>> MapVertexIndexToEdges_;

        // norms of all the edges in the mesh
        std::vector<Real3> EdgeNorms_;

        // the factor to which all the vertices on the mesh will be multiplied by
        Real factor_=1.0;
};


namespace MeshTools
{
    using Real3 = CommonTypes::Real3;
    using Real  = CommonTypes::Real;
    using index3= CommonTypes::index3;

    bool readPLY(std::string& filename, Mesh& mesh_);
    bool readPLYlibr(std::string& filename, Mesh& mesh_);
};