#include "Mesh2d.h"

namespace MeshGenerationFactory
{
    registry_<Mesh2d> registerMesh2d("Mesh2d");
}

Mesh2d::Mesh2d(ParameterPack& pack)
: MeshGeneration(pack)
{
    pack_.ReadString("file", ParameterPack::KeyType::Required, FileName_);
    pack_.ReadNumber("size_bound", ParameterPack::KeyType::Optional, size_bound_);
    pack_.ReadNumber("aspect_bound", ParameterPack::KeyType::Optional, aspect_bound_);
    isPBC_ = pack_.ReadArrayNumber("boxLength", ParameterPack::KeyType::Optional, box_length_);

    InputFileReader();

    auto& out = mesh_.accessOutput();

    for (auto name : out.getOutputNames())
    {
        output_.registerOutputFunc(name, out.getOutputFuncByName(name));
    }
}

void Mesh2d::InputFileReader()
{
    // clear the objects
    points_.clear();
    constraintIndices_.clear();

    // start reading the file
    std::ifstream ifs;
    std::stringstream ss;
    ifs.open(FileName_);

    std::string sentence;

    // Read number of points
    std::vector<Real> val;
    Real tempval;
    std::getline(ifs, sentence);
    ss.clear();
    ss.str(sentence);
    while (ss >> tempval)
    {
        val.push_back(tempval);
    }

    ASSERT((val.size() == 1), "There should be only 1 value on the first line which corresponds to number of points.");
    numPoints_ = val[0];

    // read num edges
    val.clear();
    std::getline(ifs, sentence);
    ss.clear();
    ss.str(sentence);
    while (ss >> tempval)
    {
        val.push_back(tempval);
    }

    ASSERT((val.size() == 1), "There should be only 1 value on the first line which corresponds to number of points.");
    numEdges_ = val[0];

    // read num Seeds
    val.clear();
    std::getline(ifs, sentence);
    ss.clear();
    ss.str(sentence);
    while (ss >> tempval)
    {
        val.push_back(tempval);
    }

    ASSERT((val.size() == 1), "There should be only 1 value on the first line which corresponds to number of points.");
    numSeed_ = val[0];

    // Now we have numbers, let's read points
    points_.resize(numPoints_);
    for (int i=0;i<numPoints_;i++)
    {
        std::getline(ifs, sentence);
        ss.clear();
        ss.str(sentence);
        
        std::vector<Real> tempP;
        Real2 p;
        Real val;

        while (ss >> val)
        {
            tempP.push_back(val);
        }

        ASSERT((tempP.size() == 2), "The input points should be 2 dimensional.");

        for (int j=0;j<2;j++)
        {
            p[j] = tempP[j];
        }

        points_[i] = p;
    }

    // record the number of input points
    numInputPoints_ = points_.size();

    // Now let's read edges
    constraintIndices_.resize(numEdges_);
    for (int i=0;i<numEdges_;i++)
    {
        std::getline(ifs, sentence);
        ss.clear();
        ss.str(sentence);

        std::vector<int> tempI;
        index2 idx;
        int ind;

        while (ss >> ind)
        {
            tempI.push_back(ind);
        }

        ASSERT((tempI.size() == 2), "The input edges should be 2 dimensional.");

        for (int j=0;j<2;j++)
        {
            idx[j] = tempI[j];
        }

        constraintIndices_[i] = idx;
    }

    // Now let's read the seed points
    seeds_.resize(numSeed_);
    for (int i=0;i<numSeed_;i++)
    {
        std::getline(ifs, sentence);
        ss.clear();
        ss.str(sentence);

        std::vector<Real> tempP;
        Real2 p;
        Real val;

        while (ss >> val)
        {
            tempP.push_back(val);
        }

        ASSERT((tempP.size() == 2), "The input seeds should be 2 dimensional.");

        for (int j=0;j<2;j++)
        {
            p[j] = tempP[j];
        }

        seeds_[i] = p;
    }
}

void Mesh2d::findDistance(const Real2& p1, const Real2& p2, Real2& dist, Real& distsq)
{
    for (int i=0;i<2;i++)
    {
        dist[i] = p1[i] - p2[i];
    }

    if (isPBC())
    {
        for (int i=0;i<2;i++)
        {
            if (dist[i] > box_length_[i] * 0.5) dist[i] = dist[i] - box_length_[i];
            if (dist[i] < -box_length_[i] * 0.5) dist[i] = dist[i] + box_length_[i];
        }
    }

    distsq = 0.0;
    for (int i=0;i<2;i++)
    {
        distsq += (dist[i] * dist[i]);
    }
}

void Mesh2d::checkMaxDistanceBetweenEdges()
{
    Real max = -std::numeric_limits<Real>::infinity();
    for (auto e : constraintIndices_)
    {
        Real2 p1 = points_[e[0]];
        Real2 p2 = points_[e[1]];

        Real2 distance;
        Real distsq;

        findDistance(p1, p2, distance, distsq);

        Real dist = std::sqrt(distsq);

        if (dist > max)
        {
            max = dist;
        }
    }
    MaxDistanceBetweenEdges_ = max;
    std::cout << "Max distance = " << MaxDistanceBetweenEdges_ << "\n";

    if (size_bound_ < MaxDistanceBetweenEdges_)
    {
        std::cout << "warning : the size bound is smaller than the max distance between edges with size bound = " << size_bound_ << \
        " and max distance = " << MaxDistanceBetweenEdges_ << ", this might cause boundary points to increase.\n";
    }
}

void Mesh2d::updateMesh()
{   
    auto& vertices = mesh_.accessvertices();
    auto& triangles= mesh_.accesstriangles();

    // update the vertices 
    vertices.clear();
    vertices.resize(numOutputPoints_);
    for (int i=0;i<numOutputPoints_;i++)
    {
        vertices[i].position_ = OutputVertices_[i];
    }

    // update the triangles 
    triangles.resize(numOutputFaces_);
    for (int i=0;i<numOutputFaces_;i++)
    {
        triangles[i].triangleindices_ = OutputFaces_[i];
    }

    mesh_.update();
    if (isPBC())
    {
        Real3 box2d = {{0.0,box_length_[0], box_length_[1]}};
        mesh_.setBoxLength(box2d);
    }
    mesh_.CalcVertexNormals();
}

void Mesh2d::generate()
{
    Vertex_handles_.clear();

    checkMaxDistanceBetweenEdges();

    // insert all the points
    for (auto p: points_)
    {
        Vertex_handle v = cdt_.insert(Point(p[0], p[1]));
        Vertex_handles_.push_back(v);
    }

    // insert all the vertex handles
    for (auto e : constraintIndices_)
    {
        cdt_.insert_constraint(Vertex_handles_[e[0]], Vertex_handles_[e[1]]);
    }

    // insert all the seeds 
    for (auto s : seeds_)
    {
        ListSeeds_.push_back(Point(s[0], s[1]));
    }

    CGAL::refine_Delaunay_mesh_2(cdt_, ListSeeds_.begin(), ListSeeds_.end(), Criteria(aspect_bound_, size_bound_));

    numOutputPoints_ = cdt_.number_of_vertices();
    numOutputFaces_   = 0;

    int index=0;
    for (auto it = cdt_.finite_vertices_begin(); it != cdt_.finite_vertices_end();it ++)
    {
        MapFromVertexToIndex_.insert(std::make_pair(it, index));
        Real3 point = {{0.0,it -> point()[0], it -> point()[1]}};
        OutputVertices_.push_back(point);
        index++;
    }
    // CGAL::IO::write_VTU(std::cout, cdt_, CGAL::IO::ASCII);

    OutputFaces_.clear();
    for (CDT::Finite_faces_iterator fit = cdt_.finite_faces_begin(); fit != cdt_.finite_faces_end(); fit ++)
    {
        if (fit -> is_in_domain())
        {
            index3 idx;
            idx[0] = MapFromVertexToIndex_[fit -> vertex(0)];
            idx[1] = MapFromVertexToIndex_[fit -> vertex(2)];
            idx[2] = MapFromVertexToIndex_[fit -> vertex(1)];

            OutputFaces_.push_back(idx);
            numOutputFaces_ ++;
        }
    }

    // update the mesh object 
    updateMesh();

    if (isPBC())
    {
        MakePeriodic();
    }
}

void Mesh2d::MakePeriodic()
{
    Real tol=1e-8;

    // first let's find all the boundary vertices 
    mesh_.findBoundaryVertices();
    const auto& vertices = mesh_.getvertices();
    std::vector<int> BoundaryIndices_;

    // The new vertices 
    std::vector<vertex> NewVerts;
    std::vector<int> MapOldIndexToNew(vertices.size(),-1);
    int indexTot=0;
    for (int i=0;i<vertices.size();i++)
    {
        if (mesh_.isBoundary(i))
        {
            BoundaryIndices_.push_back(i);
        }
        else
        {
            NewVerts.push_back(vertices[i]);
            MapOldIndexToNew[i] = indexTot;
            indexTot++;
        }
    }

    // ASSERT((BoundaryIndices_.size() == numInputPoints_), "The number of boundary points differs before and after \
    // CDT, before = " << numInputPoints_ << " after = " << BoundaryIndices_.size());

    // we want to check the pair-wise distances between boundary points to see if they overlap
    for (int i=0;i<BoundaryIndices_.size();i++)
    {
        int idx1 = BoundaryIndices_[i];
        std::vector<int> OverLapIndices;
        for (int j=i+1;j<BoundaryIndices_.size();j++)
        {
            int idx2 = BoundaryIndices_[j];
            Real3 distance;
            Real dist;

            mesh_.getVertexDistance(vertices[idx1], vertices[idx2], distance, dist);

            // std::cout << "Distance between " << vertices[idx1].position_[0] << " " << vertices[idx1].position_[1] << " " << vertices[idx1].position_[2] << " and " \
            // << vertices[idx2].position_[0] << " " << vertices[idx2].position_[1] << " " << vertices[idx2].position_[2] << " is " << dist << "\n";

            // if these 2 points are the same 
            if (dist <= tol)
            {
                OverLapIndices.push_back(idx2);
            }
        }

        // if there's no overlap indices 
        if (OverLapIndices.size() == 0)
        {
            if (MapOldIndexToNew[idx1] == -1)
            {
                NewVerts.push_back(vertices[idx1]);
                MapOldIndexToNew[idx1] = indexTot;
                indexTot++;
            }
        }

        for (auto id : OverLapIndices)
        {
            // first check if this particular point has already been recorded
            // then we don't push back, but we propagte the current index down
            if (MapOldIndexToNew[idx1] != -1)
            {
                MapOldIndexToNew[id] = MapOldIndexToNew[idx1];
            }
            else
            {
                NewVerts.push_back(vertices[idx1]);
                MapOldIndexToNew[idx1] = indexTot;
                MapOldIndexToNew[id] = indexTot;
                indexTot++;
            }
        }
    }

    // let's update the vertices
    auto& verts = mesh_.accessvertices();
    verts.clear();
    verts.insert(verts.end(), NewVerts.begin(), NewVerts.end());

    // let's transform the triangle indices too
    auto& tri = mesh_.accesstriangles();
    for (auto& t : tri)
    {
        for (int i=0;i<3;i++)
        {
            t.triangleindices_[i] = MapOldIndexToNew[t.triangleindices_[i]];
        }
    }

    mesh_.update();
    mesh_.CalcVertexNormals();
}