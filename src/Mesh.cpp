#include "MarchingCubesWrapper.h"

Mesh::Mesh(const ParameterPack& pack)
{
    // set up output 
    outputs_.registerOutputFunc("stl", [this](std::string name) -> void { this -> printSTL(name);});
    outputs_.registerOutputFunc("ply", [this](std::string name) -> void { this -> printPLY(name);});
    outputs_.registerOutputFunc("triangles", [this](std::string name) -> void { this -> printTriangleIndices(name);});
    outputs_.registerOutputFunc("vertex", [this](std::string name) -> void {this -> printVertices(name);});
    outputs_.registerOutputFunc("normal", [this](std::string name) -> void { this -> printNormals(name);});
    outputs_.registerOutputFunc("plyAng", [this](std::string name)-> void {this -> printPLYAng(name);});

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

void Mesh::printPLYAng(std::string name)
{
    std::cout << "Printing ply Ang" << std::endl;
    std::ofstream ofs;
    ofs.open(name);

    ofs << "ply" << "\n";
    ofs << "format ascii 1.0\n";
    ofs << "comment Created by Yusheng Cai\n";

    int sizeVertex = vertices_.size();
    int sizetriangle = triangles_.size();

    ofs << "element vertex " << sizeVertex << std::endl; 
    ofs << "property float x\n";
    ofs << "property float y\n";
    ofs << "property float z\n";
    ofs << "property float nx\n";
    ofs << "property float ny\n";
    ofs << "property float nz\n";
    ofs << "element face " << sizetriangle << "\n";
    ofs << "property list uchar uint vertex_indices\n";
    ofs << "end_header\n";

    ofs << std::fixed << std::setprecision(6);
    for (int i=0;i<sizeVertex;i++)
    {
        for (int j=0;j<3;j++)
        {
            ofs << vertices_[i].position_[j]*10.0 << " ";
        }

        for (int j=0;j<3;j++)
        {
            ofs << vertices_[i].normals_[j] << " ";
        }

        ofs << "\n";
    }

    for (int i=0;i<sizetriangle;i++)
    {
        ofs << 3 << " ";
        auto& t = triangles_[i];

        for (int j=0;j<3;j++)
        {
            ofs << t.triangleindices_[j] << " ";
        }
        ofs << "\n";
    }

    ofs.close();

}

void Mesh::printNormals(std::string name)
{
    std::ofstream ofs;
    ofs.open(name);

    for (int i=0;i<vertices_.size();i++)
    {
        for (int j=0;j<3;j++)
        {
            ofs << vertices_[i].normals_[j] << "\t";
        }
        ofs << "\n";
    }

    ofs.close();
}

void Mesh::print()
{
    for (int i=0;i<outs_.size();i++)
    {
        outputs_.getOutputFuncByName(outs_[i])(outputNames_[i]);
    }
}

void Mesh::printPLY(std::string name)
{
    std::cout << "Printing ply" << std::endl;
    std::ofstream ofs;
    ofs.open(name);

    ofs << "ply" << "\n";
    ofs << "format ascii 1.0\n";
    ofs << "comment Created by Yusheng Cai\n";

    int sizeVertex = vertices_.size();
    int sizetriangle = triangles_.size();

    ofs << "element vertex " << sizeVertex << std::endl; 
    ofs << "property float x\n";
    ofs << "property float y\n";
    ofs << "property float z\n";
    ofs << "property float nx\n";
    ofs << "property float ny\n";
    ofs << "property float nz\n";
    ofs << "element face " << sizetriangle << "\n";
    ofs << "property list uchar uint vertex_indices\n";
    ofs << "end_header\n";

    ofs << std::fixed << std::setprecision(6);
    for (int i=0;i<sizeVertex;i++)
    {
        for (int j=0;j<3;j++)
        {
            ofs << vertices_[i].position_[j] << " ";
        }

        for (int j=0;j<3;j++)
        {
            ofs << vertices_[i].normals_[j] << " ";
        }

        ofs << "\n";
    }

    for (int i=0;i<sizetriangle;i++)
    {
        ofs << 3 << " ";
        auto& t = triangles_[i];

        for (int j=0;j<3;j++)
        {
            ofs << t.triangleindices_[j] << " ";
        }
        ofs << "\n";
    }

    ofs.close();
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

        #ifdef MY_DEBUG
        std::cout << "Face " << i << " = " << index0 << " " << index1 << " " << index2 << std::endl;
        #endif

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
        #ifdef MY_DEBUG
        std::cout << "The pervetex dir 1 before anything = " << PerVertexdir1_[i][0] << " " << PerVertexdir1_[i][1] << " " << PerVertexdir1_[i][2] << std::endl;
        #endif 
        PerVertexdir1_[i] = LinAlg3x3::CrossProduct(PerVertexdir1_[i], vertices_[i].normals_);
        LinAlg3x3::normalize(PerVertexdir1_[i]);

        Real3 B = LinAlg3x3::CrossProduct(vertices_[i].normals_, PerVertexdir1_[i]) ;
        LinAlg3x3::normalize(B);

        PerVertexdir2_[i] = B;


        #ifdef MY_DEBUG
        std::cout << "For vertex " << i << " the vertex dir 1 = " << PerVertexdir1_[i][0] << " " << PerVertexdir1_[i][1] << " " << PerVertexdir1_[i][2] << "\n";
        std::cout << "For vertex " << i << " the vertex dir 2 = " << PerVertexdir2_[i][0] << " " << PerVertexdir2_[i][1] << " " << PerVertexdir2_[i][2] << "\n";
        #endif
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

// Compute per-vertex point areas
void Mesh::test()
{
    triangleArea_.clear();
    triangleArea_.resize(vertices_.size(),0.0);

    int nf = triangles_.size();

    cornerArea_.clear();
    cornerArea_.resize(nf);

	for (int i = 0; i < nf; i++) {
		// Edges
        auto& t = triangles_[i];
        Real3 edge1, edge2, edge3;
        for (int j=0;j<3;j++)
        {
            edge1[j] = vertices_[t.triangleindices_[2]].position_[j] - vertices_[t.triangleindices_[1]].position_[j];
            edge2[j] = vertices_[t.triangleindices_[0]].position_[j] - vertices_[t.triangleindices_[2]].position_[j];
            edge3[j] = vertices_[t.triangleindices_[1]].position_[j] - vertices_[t.triangleindices_[0]].position_[j];
        }

		// Compute corner weights
        Real3 cross = LinAlg3x3::CrossProduct(edge1, edge2);
        Real  area  = 0.5 * LinAlg3x3::norm(cross);

        Real e1 = LinAlg3x3::norm(edge1);
        Real e2 = LinAlg3x3::norm(edge2);
        Real e3 = LinAlg3x3::norm(edge3);

		Real3 l2 = { e1*e1, e2*e2, e3*e3};

		// Barycentric weights of circumcenter
		Real3 bcw = { l2[0] * (l2[1] + l2[2] - l2[0]),
		                 l2[1] * (l2[2] + l2[0] - l2[1]),
		                 l2[2] * (l2[0] + l2[1] - l2[2]) };

		if (bcw[0] <= 0.0f) {
			cornerArea_[i][1] = -0.25f * l2[2] * area/LinAlg3x3::DotProduct(edge1, edge3);
			cornerArea_[i][2] = -0.25f * l2[1] * area/LinAlg3x3::DotProduct(edge1, edge2);
			cornerArea_[i][0] = area - cornerArea_[i][1] - cornerArea_[i][2];
		} else if (bcw[1] <= 0.0f) {
			cornerArea_[i][2] = -0.25f * l2[0] * area/LinAlg3x3::DotProduct(edge2, edge1);
			cornerArea_[i][0] = -0.25f * l2[2] * area/LinAlg3x3::DotProduct(edge2, edge3);
			cornerArea_[i][1] = area - cornerArea_[i][2] - cornerArea_[i][0];
		} else if (bcw[2] <= 0.0f) {
            std::cout << "t " << i << " is in 3" << std::endl;
			cornerArea_[i][0] = -0.25f * l2[1] * area/LinAlg3x3::DotProduct(edge3, edge2);
			cornerArea_[i][1] = -0.25f * l2[0] * area/LinAlg3x3::DotProduct(edge3, edge1);
			cornerArea_[i][2] = area - cornerArea_[i][0] - cornerArea_[i][1];
		} else {
			float scale = 0.5f * area / (bcw[0] + bcw[1] + bcw[2]);
			for (int j = 0; j < 3; j++)
            {
                int next = j - 1;
                int nextnext = j -2;

                if (next < 0)
                {
                    next += 3;
                }

                if (nextnext < 0)
                {
                    nextnext += 3;
                }

				cornerArea_[i][j] = scale * (bcw[next] +
				                             bcw[nextnext]);
            }
		}

        #ifdef MY_DEBUG
        std::cout << "corner Area " << i << " 0 = "  << cornerArea_[i][0] << std::endl;
        std::cout << "corner Area " << i << " 1 = " << cornerArea_[i][1] << std::endl;
        std::cout << "corner Area " << i << " 2 = " << cornerArea_[i][2] << std::endl;
        #endif
		triangleArea_[t.triangleindices_[0]] += cornerArea_[i][0];
		triangleArea_[t.triangleindices_[1]] += cornerArea_[i][1];
		triangleArea_[t.triangleindices_[2]] += cornerArea_[i][2];
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

void Mesh::printVertices(std::string name)
{
    std::ofstream ofs_;

    ofs_.open(name);

    for (int i=0;i<vertices_.size();i++)
    {
        for (int j=0;j<3;j++)
        {
            ofs_ << vertices_[i].position_[j] << "\t";
        }
        ofs_ << "\n";
    }

    ofs_.close();

}

void Mesh::printTriangleIndices(std::string name)
{
    std::ofstream ofs_;

    ofs_.open(name);

    for (int i=0;i<triangles_.size();i++)
    {
        for (int j=0;j<3;j++)
        {
            ofs_ << triangles_[i].triangleindices_[j] << "\t";
        }

        ofs_ <<"\n";
    }

    ofs_.close();
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