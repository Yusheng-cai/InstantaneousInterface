#include "MarchingCubesWrapper.h"

Mesh::Mesh(const ParameterPack& pack)
{
    // set up output 
    outputs_.registerOutputFunc("stl", [this](std::string name) -> void { this -> printSTL(name);});

    // the mesh pack is the pack for the density field
    auto MeshPack = pack.findParamPack("Mesh", ParameterPack::KeyType::Optional);

    if (MeshPack != nullptr)
    {
        MeshRefineStrategyInput input = {*this, const_cast<ParameterPack&>(*MeshPack)};
        auto read = MeshPack->ReadString("type", ParameterPack::KeyType::Optional, refineStrategy_);

        if (read)
        {
            MeshRefine_ = refinePtr(MeshRefineStrategyFactory::factory::instance().create(refineStrategy_, input));
        }

        MeshPack->ReadVectorString("outputs", ParameterPack::KeyType::Optional, outs_);
        MeshPack->ReadVectorString("outputNames", ParameterPack::KeyType::Optional, outputNames_);
    }
}

void Mesh::print()
{
    for (int i=0;i<outs_.size();i++)
    {
        outputs_.getOutputFuncByName(outs_[i])(outputNames_[i]);
    }
}

void Mesh::printSTL(std::string name)
{
    std::ofstream ofs;
    ofs.open(name);

    ofs << "solid " << name << "\n";

    for (int i=0;i<triangles_.size();i++)
    {
        ofs << "facet normal " << facetNormals_[i][0] << " " << facetNormals_[i][1] << " " << facetNormals_[i][2] << "\n";

        ofs << "\touter loop\n";
        auto& t = triangles_[i];

        for (int j=0;j<3;j++)
        {
            auto& pos = t.vertices_[j].position_;
            ofs << "vertex " << pos[0] << " " << pos[1] << " " << pos[2] << "\n";
        }

        ofs << "\tendloop\n";
        ofs << "endfacet\n";
    }
    ofs << "endsolid " << name;
}

void Mesh::CalcPerVertexDir()
{
    PerVertexdir1_.resize(vertices_.size());
    PerVertexdir2_.resize(vertices_.size());

    for (int i=0;i<triangles_.size();i++)
    {
        auto& t = triangles_[i];
 
        int index0 = t.triangleindices_[0];
        int index1 = t.triangleindices_[1];
        int index2 = t.triangleindices_[2];

        Real3 diff0;
        Real3 diff1;
        Real3 diff2;

        for (int j=0;j<3;j++)
        {
            diff0[j] = vertices_[index1].position_[j] - vertices_[index0].position_[j];
            diff1[j] = vertices_[index2].position_[j] - vertices_[index1].position_[j];
            diff2[j] = vertices_[index0].position_[j] - vertices_[index2].position_[j];
        }

        PerVertexdir1_[index0] = diff0;
        PerVertexdir1_[index1] = diff1;
        PerVertexdir1_[index2] = diff2;
    }

    for (int i=0;i<PerVertexdir1_.size();i++)
    {
        PerVertexdir1_[i] = LinAlg3x3::CrossProduct(vertices_[i].normals_,PerVertexdir1_[i]);
        LinAlg3x3::normalize(PerVertexdir1_[i]);

        Real3 B = LinAlg3x3::CrossProduct(PerVertexdir1_[i],vertices_[i].normals_) ;
        LinAlg3x3::normalize(B);

        PerVertexdir2_[i] = B;
    } 
}

void Mesh::refine()
{
    if (MeshRefine_.get() != nullptr)
    {
        std::cout << "refining" << std::endl;
        MeshRefine_ -> refine();
    }
}

void Mesh::MapEdgeToFaces()
{
    // clear whatever is in the previous edgetofaces map
    MapEdgeToFace_.clear();

    // iterate over all the triangles
    for (int i=0;i<triangles_.size();i++)
    {
        // find the triangles indices
        auto& t = triangles_[i];

        for (int j=0;j<3;j++)
        {
            auto it = MapEdgeToFace_.find(t.edges_[j]);

            #ifdef MY_DEBUG
            std::cout << "For triangle " << i << " edge " << j << ", vertex1 index = " << t.edges_[j].vertex1_.index << \
            ", vertex2 index = " << t.edges_[j].vertex2_.index << std::endl;
            #endif 

            if (it == MapEdgeToFace_.end())
            {
                std::vector<int> faceIndices;
                faceIndices.push_back(i);

                MapEdgeToFace_.insert(std::make_pair(t.edges_[j], faceIndices));
            }
            else
            {
                it -> second.push_back(i);
            }
        }
    }

    // do a check to make sure that each edge is shared by at least 1 and at most 2 faces
    for (auto it=MapEdgeToFace_.begin(); it != MapEdgeToFace_.end();it++)
    {
        ASSERT((it->second.size() == 1 || it -> second.size() ==2), "The number of faces shared by an edge needs to be either 1 or 2 \
        , however the calculated result shows " << it -> second.size());
    }
}

void Mesh::findBoundaryVertices()
{
    // first get the map from edges to faces
    MapEdgeToFaces();

    for (auto it = MapEdgeToFace_.begin(); it != MapEdgeToFace_.end(); it ++)
    {
        int size = it -> second.size();
        ASSERT(( size == 1 || size == 2), "An edge can only be shared by 1 or 2 faces.");

        // if size == 1, then all the verices in the edge is part of the boundary_vertices
        if (size == 1)
        {
            bool Notfound1 = MapBoundaryVertexToBoundaryEdges_.find(it -> first.vertex1_) == MapBoundaryVertexToBoundaryEdges_.end();
            bool Notfound2 = MapBoundaryVertexToBoundaryEdges_.find(it -> first.vertex2_) == MapBoundaryVertexToBoundaryEdges_.end();

            if (Notfound1)
            {
                std::vector<edge> e;
                e.push_back(it -> first);

                MapBoundaryVertexToBoundaryEdges_.insert(std::make_pair(it -> first.vertex1_, e));
            }
            else
            {
                auto it2 = MapBoundaryVertexToBoundaryEdges_.find(it -> first.vertex1_);

                it2 -> second.push_back(it -> first);
            }

            if (Notfound2)
            {
                std::vector<edge> e;
                e.push_back(it -> first);

                MapBoundaryVertexToBoundaryEdges_.insert(std::make_pair(it -> first.vertex2_, e));
            }
            else
            {
                auto it2 = MapBoundaryVertexToBoundaryEdges_.find(it -> first.vertex2_);

                it2 -> second.push_back(it -> first);
            }
        }
    }
}

void Mesh::findVertexNeighbors()
{
    neighborIndices_.clear();
    neighborIndices_.resize(vertices_.size());

    NumNeighbors_.clear();
    NumNeighbors_.resize(vertices_.size());

    for (int i=0;i<triangles_.size();i++)
    {
        auto& t = triangles_[i];

        for (int j=0;j<3;j++)
        {
            int index1 = t.triangleindices_[j];
            for (int k=0;k<3;k++)
            {
                if (j != k)
                {
                    int index2 = t.triangleindices_[k];

                    auto it = std::find(neighborIndices_[index1].begin(), neighborIndices_[index1].end(), index2);

                    if (it == neighborIndices_[index1].end())
                    {
                        neighborIndices_[index1].push_back(index2);
                        NumNeighbors_[index1] += 1;
                    }
                }
            }
        }
    }
}

void Mesh::findTriangleIndices()
{
    VertexTriangleIndices_.clear();
    VertexTriangleIndices_.resize(vertices_.size());

    for (int i=0;i<triangles_.size();i++)
    {
        auto& t = triangles_[i];

        for (int j=0; j<3;j++)
        {
            int index = t.triangleindices_[j];

            VertexTriangleIndices_[index].push_back(i);
        }
    }
}

void Mesh::CalcTriangleAreaAndFacetNormals()
{
    triangleArea_.clear();
    triangleArea_.resize(triangles_.size());

    facetNormals_.clear();
    facetNormals_.resize(triangles_.size());

    // calculate the area of the triangles as well as the normals of the faces of triangles
    for (int i=0;i<triangles_.size();i++)
    {
        auto& t = triangles_[i];
        int index1 = t.triangleindices_[0];
        int index2 = t.triangleindices_[1];
        int index3 = t.triangleindices_[2];

        Real3 diff1 = {};
        Real3 diff2 = {};

        for (int i=0;i<3;i++)
        {
            diff1[i] = vertices_[index2].position_[i] - vertices_[index1].position_[i];
            diff2[i] = vertices_[index1].position_[i] - vertices_[index3].position_[i];
        }

        Real3 crossProduct = LinAlg3x3::CrossProduct(diff1, diff2);
        Real norm = LinAlg3x3::norm(crossProduct);
        Real Area = norm*0.5;
        Real3 normal_;
        for (int i=0;i<3;i++)
        {
            normal_[i] = crossProduct[i]/norm;
        }

        facetNormals_[i] = normal_;
        triangleArea_[i] = Area;
    }
}

void Mesh::CalcVertexNormals()
{
    // calculate vertex normals
    vertexNormals_.resize(vertices_.size());
    Real3 zeroArr = {{0,0,0}};
    std::fill(vertexNormals_.begin(), vertexNormals_.end(), zeroArr);

    // calculate the normals of each vertices
    for (int i=0;i<triangles_.size();i++)
    {
        auto& t = triangles_[i];
        Real area = triangleArea_[i];
        Real3 normals = facetNormals_[i];

        for (int j=0;j<3;j++)
        {
            int index = t.triangleindices_[j];
            auto& vNorm = vertexNormals_[index];

            for (int k=0;k<3;k++)
            {
                vNorm[k] += area*normals[k];
            }
        }
    }


    for (int i=0;i<vertexNormals_.size();i++)
    {
        Real norm = LinAlg3x3::norm(vertexNormals_[i]);

        for (int j=0;j<3;j++)
        {
            vertexNormals_[i][j] = vertexNormals_[i][j]/norm;
        }
    }

    for(int i=0;i<vertices_.size();i++)
    {
        vertices_[i].normals_ = vertexNormals_[i];
    }
}


bool MeshTools::readPLY(std::string& filename, Mesh& mesh_)
{
    std::ifstream ifs_;
    std::stringstream ss_;

    ifs_.open(filename);
    ASSERT((ifs_.is_open()), "The file with name " << filename << " is not opened.");

    auto& vertices = mesh_.accessvertices();

    std::string sentence;
    while (std::getline(ifs_, sentence))
    {
        ss_.str(sentence);
        std::string token;
        
        ss_ >> token;

        if (token.empty())
        {
            continue;
        }

        double num;
        if (StringTools::StringToType<double>(token, num))
        {
            continue;
        }
        else
        {
            std::vector<double> NumberPerLine;
            double number;
            StringTools::StringToType<double>(token, number);
            NumberPerLine.push_back(number);

            while (ss_ >> token)
            {
                bool fail = StringTools::StringToType<double>(token, number);

                ASSERT((fail == false), "The read operation failed.");
                NumberPerLine.push_back(number);
            }

            ASSERT((NumberPerLine.size()==6), "There are " << NumberPerLine.size() << " numbers on each line while it is required to have x,y,z,nx,ny,nz only");

            vertex v;

            for (int i=0;i<3;i++)
            {
                v.position_[i] = NumberPerLine[i];
            }

            for (int i=0;i<3;i++)
            {
                v.normals_[i] = NumberPerLine[i+3];
            }

            v.index = vertices.size();
            vertices.push_back(v);
        }

        ss_.clear();
    }

    return true;
}

bool MeshTools::readPLYTriangle(std::string& filename, Mesh& mesh_)
{
    std::ifstream ifs_;
    std::stringstream ss_;

    ifs_.open(filename);
    ASSERT((ifs_.is_open()), "The file with name " << filename << " is not opened.");

    auto& vertices = mesh_.getvertices();
    auto& triangles= mesh_.accesstriangles();

    std::string sentence;
    while (std::getline(ifs_, sentence))
    {
        ss_.str(sentence);
        std::string token;
        
        ss_ >> token;

        if (token.empty())
        {
            continue;
        }

        int num;
        if (StringTools::StringToType<int>(token, num))
        {
            continue;
        }
        else
        {
            // ignore the first number
            std::vector<int> NumberPerLine;
            int number;

            while (ss_ >> token)
            {
                bool fail = StringTools::StringToType<int>(token, number);

                ASSERT((fail == false), "The read operation failed.");
                NumberPerLine.push_back(number);
            }

            ASSERT((NumberPerLine.size()==3), "There are " << NumberPerLine.size() << " numbers on each line while it is required to have Id_x, Id_y, Id_z only");

            triangle t;

            for (int i=0;i<3;i++)
            {
                t.triangleindices_[i] = NumberPerLine[i];
                t.vertices_[i] = vertices[t.triangleindices_[i]];
            }

            t.edges_[0].vertex1_ = t.vertices_[0];
            t.edges_[0].vertex2_ = t.vertices_[1];

            t.edges_[1].vertex1_ = t.vertices_[1];
            t.edges_[1].vertex2_ = t.vertices_[2];

            t.edges_[2].vertex1_ = t.vertices_[2];
            t.edges_[2].vertex2_ = t.vertices_[0];

            triangles.push_back(t);
        }

        ss_.clear();
    }

    return true;
}