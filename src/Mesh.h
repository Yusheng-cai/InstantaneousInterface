#pragma once
#include "tools/Assert.h"
#include "tools/CommonTypes.h"
#include "LinAlgTools.h"
#include "tools/InputParser.h"

#include <vector>
#include <array>
#include <cmath>
#include <algorithm>
#include <sstream>

class MarchingCubesWrapper;

struct vertex
{
    using Real = CommonTypes::Real;
    using Real3= CommonTypes::Real3;

    Real3 position_;
    Real3 normals_;
};

struct triangle
{
    using index3 = CommonTypes::index3;

    index3 triangleindices_;
};

class Mesh
{
    public:
        using Real3 = CommonTypes::Real3;
        using Real  = CommonTypes::Real;

        Mesh(){};
        ~Mesh(){};

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

        int getNumVertices() const {return vertices_.size();}
        int getNumTriangles() const {return triangles_.size();}

        // Find neighbors indices for a vertex
        void findVertexNeighbors();

        // function that finds the triangle indices that a vertex belongs to 
        void findTriangleIndices();

        // function that finds the area of all the triangules on the surface
        void CalcTriangleAreaAndFacetNormals();

        // function that calculates the normals of each of the vertex
        void CalcVertexNormals();

        // Find the 
        void CalcPerVertexDir();
    
    private:
        std::vector<vertex> vertices_;
        std::vector<triangle> triangles_;

        // finds the neighbors of a vertex as well as count its number of neighbors
        std::vector<std::vector<int>> neighborIndices_;
        std::vector<int> NumNeighbors_;

        std::vector<std::vector<int>> VertexTriangleIndices_;

        std::vector<Real> triangleArea_;
        std::vector<Real3> facetNormals_;

        std::vector<Real3> vertexNormals_;

        std::vector<Real3> PerVertexdir1_;
        std::vector<Real3> PerVertexdir2_;
};


namespace MeshTools
{
    bool readPLY(std::string& filename, Mesh& mesh_);
    bool readPLYTriangle(std::string& filename, Mesh& mesh_);
};