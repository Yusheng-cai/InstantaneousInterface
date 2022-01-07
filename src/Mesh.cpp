#include "Mesh.h"

Mesh::Mesh(const ParameterPack* pack)
{
    // set up output 
    outputs_.registerOutputFunc("stl", [this](std::string name) -> void { this -> printSTL(name);});
    outputs_.registerOutputFunc("ply", [this](std::string name) -> void { this -> printPLY(name);});
    outputs_.registerOutputFunc("plyAng", [this](std::string name)-> void {this -> printPLYAng(name);});
    outputs_.registerOutputFunc("plynm", [this](std::string name) -> void {this -> printPLYnm(name);});
    outputs_.registerOutputFunc("plylibr", [this](std::string name) -> void {this -> printPLYlibr(name);});
    outputs_.registerOutputFunc("boundary", [this](std::string name) -> void {this -> printBoundaryVertices(name);});
    outputs_.registerOutputFunc("area", [this](std::string name) -> void {this -> printArea(name);});


    if (pack != nullptr)
    {
        MeshRefineStrategyInput input = {*this, const_cast<ParameterPack&>(*pack)};
        auto read = pack->ReadString("type", ParameterPack::KeyType::Optional, refineStrategy_);

        if (read)
        {
            MeshRefine_ = refinePtr(MeshRefineStrategyFactory::factory::instance().create(refineStrategy_, input));
        }

        pack->ReadVectorString("outputs", ParameterPack::KeyType::Optional, outs_);
        pack->ReadVectorString("outputNames", ParameterPack::KeyType::Optional, outputNames_);
        isPeriodic_ = pack->ReadArrayNumber("BoxLength", ParameterPack::KeyType::Optional, boxLength_);
        pack->ReadNumber("factor", ParameterPack::KeyType::Optional, factor_);
    }
}

void Mesh::printArea(std::string name)
{
    CalcTriangleAreaAndFacetNormals();

    std::ofstream ofs;
    ofs.open(name);

    for (int i=0;i<triangleArea_.size();i++)
    {
        ofs << triangleArea_[i] << "\n";
    }

    ofs.close();
}

void Mesh::printBoundaryVertices(std::string name)
{
    std::ofstream ofs;
    ofs.open(name);

    MapEdgeToFaces();

    findBoundaryVertices();

    for (int i=0;i<vertices_.size();i++)
    {
        if (isBoundary(i))
        {
            ofs << i << "\n";
        }
    }

    ofs.close();
}

void Mesh::printPLYlibr(std::string name)
{
    std::vector<std::array<double,3>> positions(vertices_.size());
    std::vector<std::vector<size_t>> fInd(triangles_.size());
    std::vector<std::vector<double>> normals(3, std::vector<double>(vertices_.size(),0.0));

    std::vector<std::string> directioNames = {"nx", "ny", "nz"};

    for (int i=0;i<vertices_.size();i++)
    {
        for (int j=0;j<3;j++)
        {
            positions[i][j] = factor_ * vertices_[i].position_[j];
            normals[j][i] = vertices_[i].normals_[j];
        }
    }

    for (int i=0;i<triangles_.size();i++)
    {
        for (int j=0;j<3;j++)
        {
            fInd[i].push_back(triangles_[i].triangleindices_[j]);
        }
    }

    happly::PLYData plyOut;

    plyOut.addVertexPositions(positions);
    plyOut.addFaceIndices(fInd);

    for (int i=0;i<3;i++)
    {
        plyOut.getElement("vertex").addProperty(directioNames[i], normals[i]);
    }

    plyOut.write(name, happly::DataFormat::ASCII);
    std::cout << "Done reading." << std::endl;
}

void Mesh::printPLYnm(std::string name)
{
    std::cout << "Printing ply nm" << std::endl;
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
            ofs << vertices_[i].position_[j]/10.0 << " ";
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

        Real3 diff0, diff1, diff2;
        Real diffsq0, diffsq1, diffsq2;

        getVertexDistance(vertices_[index1], vertices_[index0], diff0, diffsq0);
        getVertexDistance(vertices_[index2], vertices_[index1], diff1, diffsq1);
        getVertexDistance(vertices_[index0], vertices_[index2], diff2, diffsq2);

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

        // map vertex to edges 
        for (int j=0;j<3;j++)
        {
            int index1 = t.edges_[j].vertex1_.index;
            int index2 = t.edges_[j].vertex2_.index;

            std::array<int,2> indices = {{index1,index2}};

            for (int k=0;k<2;k++)
            {
                auto it = MapVertexIndexToEdges_.find(indices[k]); 

                if (it == MapVertexIndexToEdges_.end())
                {
                    std::vector<edge> temp;
                    temp.push_back(t.edges_[j]);
                    MapVertexIndexToEdges_.insert(std::make_pair(indices[k], temp));
                }
                else
                {
                    // find if the edge is already in the vector
                    auto f = std::find(it ->second.begin(), it ->second.end(), t.edges_[j]);

                    if (f == it -> second.end())
                    {
                        it -> second.push_back(t.edges_[j]);
                    }
                }
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

void Mesh::MapVertexToFaces()
{
    MapVertexIndicesToFaceIndices_.clear();
    MapVertexIndicesToFaceIndices_.resize(vertices_.size());

    for (int i=0;i<triangles_.size();i++)
    {
        auto t = triangles_[i].triangleindices_;
        for (int j=0;j<3;j++)
        {
            MapVertexIndicesToFaceIndices_[t[j]].push_back(i);
        }
    }

}

std::vector<int>& Mesh::getFaceIndicesForEdge(const edge& e)
{
    auto it = MapEdgeToFace_.find(e);

    ASSERT((it != MapEdgeToFace_.end()), "The edge with vertex indices " << e.vertex1_.index << " and " << e.vertex2_.index << " does not exist in the map from edge to face.");

    return it -> second;
}

std::vector<edge>& Mesh::getEdgeForVertex(int i)
{
    auto it = MapVertexIndexToEdges_.find(i);

    ASSERT((it != MapVertexIndexToEdges_.end()), "The vertex index " << i << " is not found with mapping to edge");

    return it -> second;
}

void Mesh::findBoundaryVertices()
{
    for (auto it = MapEdgeToFace_.begin(); it != MapEdgeToFace_.end(); it ++)
    {
        int size = it -> second.size();
        ASSERT(( size == 1 || size == 2), "An edge can only be shared by 1 or 2 faces.");

        // if size == 1, then all the vertices in the edge is part of the boundary_vertices
        if (size == 1)
        {
            bool Notfound1 = MapBoundaryVertexToBoundaryEdges_.find(it -> first.vertex1_) == MapBoundaryVertexToBoundaryEdges_.end();
            bool Notfound2 = MapBoundaryVertexToBoundaryEdges_.find(it -> first.vertex2_) == MapBoundaryVertexToBoundaryEdges_.end();

            if (Notfound1)
            {
                std::vector<edge> e;
                e.push_back(it -> first);

                MapBoundaryVertexToBoundaryEdges_.insert(std::make_pair(it -> first.vertex1_, e));
                MapBoundaryVertexIndicesToTrue_.insert(std::make_pair(it->first.vertex1_.index,true));
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
                MapBoundaryVertexIndicesToTrue_.insert(std::make_pair(it->first.vertex2_.index,true));
            }
            else
            {
                auto it2 = MapBoundaryVertexToBoundaryEdges_.find(it -> first.vertex2_);

                it2 -> second.push_back(it -> first);
            }
        }
    }
}

bool Mesh::isBoundary(int i)
{
    auto it = MapBoundaryVertexIndicesToTrue_.find(i);

    if (it == MapBoundaryVertexIndicesToTrue_.end())
    {
        return false;
    }
    else
    {
        return true;
    }

    return false;
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

Mesh::Real Mesh::calculateVolume()
{
    Real volume_ =0.0;
    for (int i=0;i<triangles_.size();i++)
    {
        auto& t = triangles_[i];

        Real3 g;
        g.fill(0);
        for (int j=0;j<3;j++)
        {
            int index= t.triangleindices_[j];
            for (int k=0;k<3;k++)
            {
                g[k] += vertices_[index].position_[k] /3.0;
            }
        }

        Real3 vec1, vec2;
        int index1 = t.triangleindices_[0];
        int index2 = t.triangleindices_[1];
        int index3 = t.triangleindices_[2];
        for (int j=0;j<3;j++)
        {
            vec1[j] = vertices_[index2].position_[j] - vertices_[index1].position_[j];
            vec2[j] = vertices_[index3].position_[j] - vertices_[index1].position_[j];
        }

        Real3 N = LinAlg3x3::CrossProduct(vec1, vec2);

        volume_ += 1.0/6.0 * LinAlg3x3::DotProduct(g,N);
    }

    return volume_;
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

void Mesh::scaleVertices(Real num)
{
    #pragma omp parallel for
    for (int i=0;i<vertices_.size();i++)
    {
        for (int j=0;j<3;j++)
        {
            vertices_[i].position_[j] = num * vertices_[i].position_[j];
        }
    }

    update();
}
void Mesh::getVertexDistance(const vertex& v1, const vertex& v2, Real3& distVec, Real& dist)
{
    for (int i=0;i<3;i++)
        {
            Real diff = v1.position_[i] - v2.position_[i];
            distVec[i] = diff;
    }

    dist = LinAlg3x3::norm(distVec);

    if (isPeriodic())
    {
        for (int i=0;i<3;i++)
        {
            if (distVec[i] > boxLength_[i] * 0.5) distVec[i] -= boxLength_[i];
            if (distVec[i] < -boxLength_[i] * 0.5) distVec[i] += boxLength_[i];
        }

        dist = LinAlg3x3::norm(distVec);
    }
}

void Mesh::CalcTriangleAreaAndFacetNormals()
{
    triangleArea_.clear();
    triangleArea_.resize(triangles_.size());

    facetNormals_.clear();
    facetNormals_.resize(triangles_.size());

    EdgeNorms_.clear();
    EdgeNorms_.resize(triangles_.size());

    // calculate the area of the triangles as well as the normals of the faces of triangles
    for (int i=0;i<triangles_.size();i++)
    {
        auto& t = triangles_[i];
        int index1 = t.triangleindices_[0];
        int index2 = t.triangleindices_[1];
        int index3 = t.triangleindices_[2];

        Real3 diff1 = {};
        Real3 diff2 = {};
        Real3 diff3 = {};
        Real norm1, norm2, norm3;

        getVertexDistance(vertices_[index1], vertices_[index2], diff1, norm1);
        getVertexDistance(vertices_[index2], vertices_[index3], diff2, norm2);
        getVertexDistance(vertices_[index3], vertices_[index1], diff3, norm3);

        // update the edge norms of a triangle
        EdgeNorms_[i][0] = norm1;
        EdgeNorms_[i][1] = norm2;
        EdgeNorms_[i][2] = norm3;

        Real3 crossProduct = LinAlg3x3::CrossProduct(diff1, diff2);
        Real norm = LinAlg3x3::norm(crossProduct);
        Real Area = norm*0.5;

        Real3 normal_;
        for (int j=0;j<3;j++)
        {
            normal_[j] = crossProduct[j]/norm;
        }

        facetNormals_[i] = normal_;
        triangleArea_[i] = Area;
    }
}

void Mesh::updateNormals()
{
    CalcTriangleAreaAndFacetNormals();

    // clear the normals in the vertices 
    for (int i=0;i<vertices_.size();i++)
    {
        vertices_[i].normals_.fill(0);
    }

    for (int i=0;i<triangles_.size();i++)
    {
        auto& t = triangles_[i].triangleindices_;
        Real factor1 = 1.0/(EdgeNorms_[i][0] * EdgeNorms_[i][2]);
        Real factor2 = 1.0/(EdgeNorms_[i][0] * EdgeNorms_[i][1]);
        Real factor3 = 1.0/(EdgeNorms_[i][1] * EdgeNorms_[i][2]);

        int index1 = t[0];
        int index2 = t[1];
        int index3 = t[2];

        for (int j=0;j<3;j++)
        {
            vertices_[index1].normals_[j] += factor1 * facetNormals_[i][j];
            vertices_[index2].normals_[j] += factor2 * facetNormals_[i][j];
            vertices_[index3].normals_[j] += factor3 * facetNormals_[i][j];
        }
    }

    for (int i=0;i<vertices_.size();i++)
    {
        LinAlg3x3::normalize(vertices_[i].normals_);
    }

    update();
}

Mesh::Real3 Mesh::getShiftedVertexPosition(const vertex& v1, const vertex& v2)
{
    Real3 dist;
    Real distsq;

    // v1 - v2
    getVertexDistance(v1, v2, dist, distsq);
    Real3 ret;

    for (int i=0;i<3;i++)
    {
        ret[i] = v2.position_[i] + dist[i];
    }

    return ret;
}

// Compute per-vertex point areas
void Mesh::CalculateCornerArea()
{
    triangleArea_.clear();
    triangleArea_.resize(vertices_.size(),0.0);

    int nf = triangles_.size();

    cornerArea_.clear();
    cornerArea_.resize(nf);

	for (int i = 0; i < nf; i++) {
		// Edges
        auto& t = triangles_[i];
        int index1 = t.triangleindices_[0];
        int index2 = t.triangleindices_[1];
        int index3 = t.triangleindices_[2];

        Real3 edge1, edge2, edge3;
        Real edge1sq, edge2sq, edge3sq;

        getVertexDistance(vertices_[index3], vertices_[index2], edge1, edge1sq);
        getVertexDistance(vertices_[index1], vertices_[index3], edge2, edge2sq);
        getVertexDistance(vertices_[index2], vertices_[index1], edge3, edge3sq);

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

void Mesh::update()
{
    for (int i=0;i<triangles_.size();i++)
    {
        auto& t = triangles_[i];

        for (int j=0;j<3;j++)
        {
            int index = t.triangleindices_[j];
            t.vertices_[j] = vertices_[index];
        }

        // construct the edges of each triangle
        t.edges_[0].vertex1_ = t.vertices_[0];
        t.edges_[0].vertex2_ = t.vertices_[1];

        t.edges_[1].vertex1_ = t.vertices_[1];
        t.edges_[1].vertex2_ = t.vertices_[2];

        t.edges_[2].vertex1_ = t.vertices_[2];
        t.edges_[2].vertex2_ = t.vertices_[0];
    }
}

void Mesh::CalcVertexNormals()
{
    // First calculate triangle faces normals 
    CalcTriangleAreaAndFacetNormals();

    // calculate vertex normals
    vertexNormals_.resize(vertices_.size());
    Real3 zeroArr = {{0,0,0}};
    std::fill(vertexNormals_.begin(), vertexNormals_.end(), zeroArr);

    // calculate the normals of each vertices
    for (int i=0;i<triangles_.size();i++)
    {
        auto& t = triangles_[i];
        Real3 normals = facetNormals_[i];

        for (int j=0;j<3;j++)
        {
            int index = t.triangleindices_[j];
            auto& vNorm = vertexNormals_[index];

            for (int k=0;k<3;k++)
            {
                vNorm[k] += normals[k];
            }
        }
    }


    for (int i=0;i<vertexNormals_.size();i++)
    {
        Real norm = LinAlg3x3::norm(vertexNormals_[i]);

        if (norm < 1e-5)
        {
            std::cout << " vertex " << i << " " <<  vertexNormals_[i][0] << " " << vertexNormals_[i][1] << " " << vertexNormals_[i][2] << "\n";
        }

        for (int j=0;j<3;j++)
        {
            vertexNormals_[i][j] = vertexNormals_[i][j]/norm;
        }
    }

    for (int i=0;i<vertices_.size();i++)
    {
        vertices_[i].normals_ = vertexNormals_[i];
    }

}

void Mesh::CalcVertexNormalsAreaWeighted()
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

    for (int i=0;i<vertices_.size();i++)
    {
        vertices_[i].normals_ = vertexNormals_[i];
    }
}

bool MeshTools::readPLYlibr(std::string& filename, Mesh& mesh_)
{
    happly::PLYData plydata(filename);
    std::vector<std::array<double,3>> vPos = plydata.getVertexPositions();
    std::vector<std::vector<size_t>> fInd = plydata.getFaceIndices<size_t>();
    std::vector<std::string> normalNames = {"nx", "ny", "nz"};
    std::vector<bool> hasProperty_(3, false);
    std::vector<Real3> normals_;
    normals_.resize(vPos.size());

    auto names = plydata.getElement("vertex").getPropertyNames();
    bool hasnormals=false;
    for (int i=0;i<3;i++)
    {
        hasProperty_[i] = std::find(names.begin(), names.end(), normalNames[i]) != names.end();
    }

    if (hasProperty_[0] && hasProperty_[1] && hasProperty_[2])
    {
        hasnormals=true;
    }
    std::cout << "has normals = " << hasnormals << std::endl;

    // do something when it does have normals
    if (hasnormals)
    {
        for (int i=0;i<3;i++)
        {
            const auto& n = plydata.getElement("vertex").getProperty<double>(normalNames[i]);
            ASSERT((n.size() == vPos.size()), "The size of nx must equal to the number of vertices.");

            for (int j=0;j<vPos.size();j++)
            {
                normals_[j][i] = n[j];
            }
        }
    }

    // first we update the mesh
    auto& vertices = mesh_.accessvertices();
    auto& triangles= mesh_.accesstriangles();
    vertices.clear();
    triangles.clear();

    vertices.resize(vPos.size());
    triangles.resize(fInd.size());

    for (int i=0;i<vertices.size();i++)
    {
        auto& v = vertices[i];
        for (int j=0;j<3;j++)
        {
            v.position_[j] = vPos[i][j];
            v.normals_[j]  = normals_[i][j];
        }
        v.index = i;
    }

    for (int i=0;i<triangles.size();i++)
    {
        auto& t = triangles[i];

        ASSERT((fInd[i].size() == 3), "We are reading triangles here while the provided face is " << fInd[i].size());

        for (int j=0;j<3;j++)
        {
            t.triangleindices_[j] = fInd[i][j];
            t.vertices_[j] = vertices[fInd[i][j]];
        }

        t.edges_[0].vertex1_ = t.vertices_[0];
        t.edges_[0].vertex2_ = t.vertices_[1];

        t.edges_[1].vertex1_ = t.vertices_[1];
        t.edges_[1].vertex2_ = t.vertices_[2];

        t.edges_[2].vertex1_ = t.vertices_[2];
        t.edges_[2].vertex2_ = t.vertices_[0];
    }

    // calculate triangle areas and facet normals 
    mesh_.CalcTriangleAreaAndFacetNormals();

    if ( ! hasnormals)
    {
        std::cout << "Calculating normals by myself." << std::endl;
        mesh_.CalcVertexNormals();

        mesh_.update();
    }

    return true;
}

bool MeshTools::readPLY(std::string& filename, Mesh& mesh_)
{
    auto& vertices = mesh_.accessvertices();
    auto& triangles= mesh_.accesstriangles();
    vertices.clear();
    triangles.clear();

    // open the file
    std::ifstream ifs_;
    std::stringstream ss_;
    ifs_.open(filename);

    if (! ifs_.is_open())
    {
        return false;
    }

    std::string sentence;
    int numfaces;
    int numvertex;
    while (std::getline(ifs_, sentence))
    {
        ss_.str(sentence);
        std::string token;

        std::vector<std::string> vectorstring;
        while (ss_ >> token)
        {
            vectorstring.push_back(token);
        }

        if (vectorstring[0] == "end_header")
        {
            break;
        }

        if (vectorstring[0] == "format")
        {
            ASSERT((vectorstring[1] == "ascii"), "Currently do not support non-ascii format ply files.");
        }

        if (vectorstring[0] == "element")
        {
            ASSERT((vectorstring.size() == 3), "The element can only have 3.");
            if (vectorstring[1] == "vertex")
            {
                numvertex = StringTools::StringToType<int>(vectorstring[2]);
            }

            if (vectorstring[1] == "face")
            {
                numfaces = StringTools::StringToType<int>(vectorstring[2]);
            }
        }

        ss_.clear();
    }

    ss_.clear();
    // read in the vertices as well as their normals     
    std::string datasentence;
    for (int i=0;i<numvertex;i++)
    {
        std::getline(ifs_, datasentence);
        ss_.str(datasentence);
        std::string data;
        std::vector<std::string> vectordata;

        vertex v;

        while (ss_ >> data)
        {
            vectordata.push_back(data);
        }

        ASSERT((vectordata.size() == 6), "The ply file must contain x, y, z, nx, ny ,nz");

        Real3 position;
        Real3 normals;

        for (int j=0;j<3;j++)
        {
            position[j] = StringTools::StringToType<Real>(vectordata[j]);
            normals[j] = StringTools::StringToType<Real>(vectordata[j+3]);
        }

        v.position_ = position;
        v.normals_ = normals;
        v.index = i;

        vertices.push_back(v);

        ss_.clear();
    }

    ss_.clear();
    // read in the triangles
    std::string trianglesentence_;
    for (int i=0;i<numfaces;i++)
    {
        std::getline(ifs_,trianglesentence_);
        ss_.str(trianglesentence_);
        std::string data;
        std::vector<std::string> vectordata;

        while (ss_ >> data)
        {
            vectordata.push_back(data);
        }

        ASSERT((vectordata.size() == 4), "The triangle file must contain 3 f1 f2 f3");

        triangle t;

        index3 faceid;

        for (int j=0;j<3;j++)
        {
            faceid[j] = StringTools::StringToType<int>(vectordata[j+1]);
        }

        t.triangleindices_ = faceid;
        for (int j=0;j<3;j++)
        {
            t.vertices_[j] = vertices[faceid[j]];
        }

        t.edges_[0].vertex1_ = t.vertices_[0];
        t.edges_[0].vertex2_ = t.vertices_[1];

        t.edges_[1].vertex1_ = t.vertices_[1];
        t.edges_[1].vertex2_ = t.vertices_[2];

        t.edges_[2].vertex1_ = t.vertices_[2];
        t.edges_[2].vertex2_ = t.vertices_[0];

        triangles.push_back(t);
        ss_.clear();
    }

    ifs_.close();

    return true;
}