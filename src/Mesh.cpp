#include "Mesh.h"

Mesh::Mesh(const ParameterPack* pack)
: pack_(pack)
{
    registerFunc();

    if (pack != nullptr)
    {
        pack->ReadVectorString("outputs", ParameterPack::KeyType::Optional, outs_);
        pack->ReadVectorString("outputNames", ParameterPack::KeyType::Optional, outputNames_);
        pack->ReadNumber("factor", ParameterPack::KeyType::Optional, factor_);

        // if box length is read, then the mesh is periodic 
        isPeriodic_ = pack->ReadArrayNumber("BoxLength", ParameterPack::KeyType::Optional, boxLength_);
    }
}

void Mesh::registerFunc()
{   
    // set up output 
    outputs_.registerOutputFunc("stl", [this](std::string name) -> void { this -> printSTL(name);});
    // print ply file by myself 
    outputs_.registerOutputFunc("ply", [this](std::string name) -> void { this -> printPLY(name);});
    // print ply file but using library
    outputs_.registerOutputFunc("boundary", [this](std::string name) -> void {this -> printBoundaryVertices(name);});
    outputs_.registerOutputFunc("area", [this](std::string name) -> void {this -> printArea(name);});
    outputs_.registerOutputFunc("cutted", [this](std::string name) -> void {this -> printCuttedMesh(name);});
    outputs_.registerOutputFunc("nonpbcMesh", [this](std::string name) -> void {this -> printNonPBCMesh(name);});
    outputs_.registerOutputFunc("translate",[this](std::string name) -> void {this -> printTranslatedMesh(name);});
    outputs_.registerOutputFunc("neighbor", [this](std::string name) -> void {this -> printNeighbors(name);});
    outputs_.registerOutputFunc("NonPeriodicTriangle", [this](std::string name) -> void {this -> printNonPeriodicTriangleIndices(name);});
}

void Mesh::printNeighbors(std::string name)
{
    // calculate the neighbor indices 
    std::vector<std::vector<int>> neighborInd;
    MeshTools::CalculateVertexNeighbors(*this, neighborInd);

    std::ofstream ofs;
    ofs.open(name);

    for (int i=0;i<neighborInd.size();i++)
    {
        ofs << i << " ";
        for (auto ind : neighborInd[i])
        {
            ofs << ind << " ";
        }
        ofs << "\n";
    }

    ofs.close();
}

void Mesh::printArea(std::string name)
{
    // calculate the area and facet normals 
    std::vector<Real> triangleA;
    std::vector<Real3> facetN;
    MeshTools::CalculateTriangleAreasAndFaceNormals(*this, triangleA, facetN);

    std::ofstream ofs;
    ofs.open(name);

    for (int i=0;i<triangleA.size();i++)
    {
        ofs << triangleA[i] << "\n";
    }

    ofs.close();
}

void Mesh::printNonPeriodicTriangleIndices(std::string name)
{
    std::ofstream ofs;
    ofs.open(name);
    std::vector<int> NonPeriodicTriangleIndices;

    if (isPeriodic())
    {
        for (int i =0;i<triangles_.size();i++)
        {
            auto t = triangles_[i].triangleindices_;
            if (MeshTools::IsPeriodicTriangle(vertices_, t, getBoxLength()))
            {
                NonPeriodicTriangleIndices.push_back(i);
            }
        }
    }
    else
    {
        NonPeriodicTriangleIndices.resize(triangles_.size());
        std::iota(NonPeriodicTriangleIndices.begin(), NonPeriodicTriangleIndices.end(),0);
    }

    for (int index : NonPeriodicTriangleIndices)
    {
        ofs << index << " ";
    }

    ofs.close();
}

void Mesh::printTranslatedMesh(std::string name)
{
    Real3 translate;
    bool readTranslate = pack_->ReadArrayNumber("translation", ParameterPack::KeyType::Optional, translate);

    if (readTranslate)
    {
        std::vector<Real3> tempVertices;
        std::vector<INT3> tempFaces;
        tempVertices.resize(vertices_.size());
        int index=0;
        for (auto& v : vertices_)
        {
            for (int i=0;i<3;i++)
            {
                 tempVertices[index][i] = v.position_[i] + translate[i];
            }
            index++;
        }

        for (auto& t : triangles_)
        {
            tempFaces.push_back(t.triangleindices_);
        }
        MeshTools::writePLY(name, tempVertices, tempFaces, factor_);
    }
}

void Mesh::printNonPBCMesh(std::string name)
{
    if (isPeriodic())
    {
        std::vector<Real3> tempVertices;
        std::vector<INT3> tempFaces;
        MeshTools::ConvertToNonPBCMesh(*this, tempVertices, tempFaces);

        MeshTools::writePLY(name, tempVertices, tempFaces, factor_);
    }
    else
    {
        std::cout << "WARNING: You asked to print NON PBC Mesh while the Mesh is not PBC." << "\n";
    }
}

void Mesh::MoveVertexIntoBox(const Real3& OldVerPos, Real3& NewVerPos)
{
    if (isPeriodic())
    {
        Real3 boxCenter;
        for (int i=0;i<3;i++)
        {
            boxCenter[i] = boxLength_[i] * 0.5;
        }

        Real3 diff;
        for (int i=0;i<3;i++)
        {
            diff[i] = OldVerPos[i] - boxCenter[i];

            if (diff[i] >= boxLength_[i] * 0.5) diff[i] = diff[i] - boxLength_[i];
            if (diff[i] <= -boxLength_[i] * 0.5) diff[i] = diff[i] + boxLength_[i]; 

            NewVerPos[i] = boxCenter[i] + diff[i];
        }
    }
}

void Mesh::printBoundaryVertices(std::string name)
{
    std::ofstream ofs;
    ofs.open(name);

    std::vector<bool> boundaryIndicator;
    MeshTools::MapEdgeToFace(*this, MapEdgeToFace_, MapVertexToEdges_);
    MeshTools::CalculateBoundaryVertices(*this, MapEdgeToFace_, boundaryIndicator);

    ofs << "# Index Vx Vy Vz Nx Ny Nz" << "\n";
    for (int i=0;i<vertices_.size();i++)
    {
        if (MeshTools::IsBoundary(i, boundaryIndicator))
        {
            ofs << i << " " << vertices_[i].position_[0] << " " << vertices_[i].position_[1] << " " << vertices_[i].position_[2] << \
            " " << vertices_[i].normals_[0] << " " << vertices_[i].normals_[1] << " " << vertices_[i].normals_[2] << "\n";
        }
    }

    ofs.close();
}

void Mesh::printCuttedMesh(std::string name)
{
    ASSERT((pack_ != nullptr), "If you wanted to print cutted Mesh, you must provide a parameter pack for mesh.");
    Real3 volume;

    pack_->ReadArrayNumber("largerthan", ParameterPack::KeyType::Required, volume);

    int index=0;
    std::map<int,int> MapOldIndexToNew;
    std::vector<Real3> newVertices;
    for (int i=0;i<vertices_.size();i++)
    {
        auto& v = vertices_[i];
        if (v.position_[0] >= volume[0] && v.position_[1] >= volume[1] && v.position_[2] >= volume[2])
        {
            int newindex = index;
            newVertices.push_back(v.position_);
            MapOldIndexToNew.insert(std::make_pair(i, newindex));
            index ++;
        }
    }

    std::vector<INT3> newTriangles;
    for (auto& t : triangles_)
    {
        bool it1 = MapOldIndexToNew.find(t.triangleindices_[0]) != MapOldIndexToNew.end();
        bool it2 = MapOldIndexToNew.find(t.triangleindices_[1]) != MapOldIndexToNew.end();
        bool it3 = MapOldIndexToNew.find(t.triangleindices_[2]) != MapOldIndexToNew.end();

        if (it1 && it2 && it3)
        {
            INT3 newt;
            for (int i=0;i<3;i++)
            {
                auto it = MapOldIndexToNew.find(t.triangleindices_[i]);
                int newIndex = it -> second;
                newt[i] = newIndex;
            }
            newTriangles.push_back(newt);
        }
    }
    MeshTools::writePLY(name, newVertices, newTriangles);
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
            ofs << vertices_[i].position_[j] * factor_ << " ";
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
        PerVertexdir1_[i] = LinAlg3x3::CrossProduct(PerVertexdir1_[i], vertices_[i].normals_);
        LinAlg3x3::normalize(PerVertexdir1_[i]);

        Real3 B = LinAlg3x3::CrossProduct(vertices_[i].normals_, PerVertexdir1_[i]) ;
        LinAlg3x3::normalize(B);

        PerVertexdir2_[i] = B;
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

void Mesh::getVertexDistance(const Real3& v1, const Real3& v2, Real3& distVec, Real& dist)
{
    for (int i=0;i<3;i++)
        {
            Real diff = v1[i] - v2[i];
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

void Mesh::getVertexDistance(const vertex& v1, const vertex& v2, Real3& distVec, Real& dist)
{
    getVertexDistance(v1.position_, v2.position_, distVec, dist);
}

void Mesh::CalculateShift(const Real3& v1, const Real3& v2, Real3& shiftVec)
{
    for (int i=0;i<3;i++)
    {
        Real d = v1[i] - v2[i];

        if (d > (0.5 * boxLength_[i])) shiftVec[i] = -boxLength_[i];
        else if (d < (- 0.5 * boxLength_[i])) shiftVec[i] = boxLength_[i];
        else shiftVec[i] = 0.0;
    }
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

Mesh::Real3 Mesh::getShiftIntoBox(const Real3& v1)
{
    Real3 shift;
    Real3 center;
    for (int i=0;i<3;i++)
    {
        center[i] = boxLength_[i] * 0.5;
    }

    for (int i=0;i<3;i++)
    {
        Real diff = v1[i] - center[i];

        if (diff >= boxLength_[i] * 0.5) shift[i] = -boxLength_[i];
        else if (diff <= - boxLength_[i] * 0.5) shift[i] = boxLength_[i];
        else shift[i] = 0.0;
    }

    return shift;
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
    }
}

void Mesh::CalcVertexNormals()
{
    MeshTools::CalculateTriangleAreasAndFaceNormals(*this, triangleArea_, facetNormals_);

    // calculate vertex normals
    std::vector<Real3> VertNorms(getNumVertices(), {{0,0,0}});

    // calculate the normals of each vertices
    for (int i=0;i<triangles_.size();i++)
    {
        auto& t = triangles_[i];
        Real3 normals = facetNormals_[i];

        for (int j=0;j<3;j++)
        {
            int index = t.triangleindices_[j];

            for (int k=0;k<3;k++)
            {
                VertNorms[index][k] += normals[k];
            }
        }
    }

    #pragma omp parallel for
    for (int i=0;i<getNumVertices();i++)
    {
        Real norm = LinAlg3x3::norm(VertNorms[i]);

        for (int j=0;j<3;j++)
        {
            vertices_[i].normals_[j] = VertNorms[i][j]/norm;
        }
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

                                    /**************************************
                                     ************   MeshTools *************
                                     *************************************/       

void MeshTools::writePLY(std::string filename, const std::vector<Real3>& vertices, const std::vector<INT3>& faces, const std::vector<Real3>& normals)
{
    std::ofstream ofs;
    ofs.open(filename);

    ofs << "ply" << "\n";
    ofs << "format ascii 1.0\n";
    ofs << "comment Created by Yusheng Cai\n";

    int sizeVertex = vertices.size();
    int sizetriangle = faces.size();

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
            ofs << vertices[i][j] << " ";
        }

        for (int j=0;j<3;j++)
        {
            ofs << normals[i][j] << " ";
        }

        ofs << "\n";
    }

    for (int i=0;i<sizetriangle;i++)
    {
        ofs << 3 << " ";
        auto& t = faces[i];

        for (int j=0;j<3;j++)
        {
            ofs << t[j] << " ";
        }
        ofs << "\n";
    }

    ofs.close();
}


void MeshTools::writePLY(std::string filename, const std::vector<Real3>& vertices, const std::vector<INT3>& faces, \
Real factor)
{
    std::ofstream ofs;
    ofs.open(filename);

    ofs << "ply" << "\n";
    ofs << "format ascii 1.0\n";
    ofs << "comment Created by Yusheng Cai\n";

    int sizeVertex = vertices.size();
    int sizetriangle = faces.size();

    ofs << "element vertex " << sizeVertex << std::endl; 
    ofs << "property float x\n";
    ofs << "property float y\n";
    ofs << "property float z\n";
    ofs << "element face " << sizetriangle << "\n";
    ofs << "property list uchar uint vertex_indices\n";
    ofs << "end_header\n";

    ofs << std::fixed << std::setprecision(6);
    for (int i=0;i<sizeVertex;i++)
    {
        for (int j=0;j<3;j++)
        {
            ofs << vertices[i][j] * factor << " ";
        }

        ofs << "\n";
    }

    for (int i=0;i<sizetriangle;i++)
    {
        ofs << 3 << " ";
        auto& t = faces[i];

        for (int j=0;j<3;j++)
        {
            ofs << t[j] << " ";
        }
        ofs << "\n";
    }

    ofs.close();
}

bool MeshTools::readPLYlibr(std::string& filename, Mesh& mesh)
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
    auto& vertices = mesh.accessvertices();
    auto& triangles= mesh.accesstriangles();
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
    }

    if ( ! hasnormals)
    {
        std::cout << "Calculating normals by myself." << std::endl;
        mesh.CalcVertexNormals();

        mesh.update();
    }

    return true;
}

bool MeshTools::readPLY(std::string& filename, Mesh& mesh)
{
    auto& vertices = mesh.accessvertices();
    auto& triangles= mesh.accesstriangles();
    vertices.clear();
    triangles.clear();

    // open the file
    std::ifstream ifs;
    std::stringstream ss;
    ifs.open(filename);

    if (! ifs.is_open())
    {
        return false;
    }

    std::string sentence;
    int numfaces;
    int numvertex;
    while (std::getline(ifs, sentence))
    {
        ss.str(sentence);
        std::string token;

        std::vector<std::string> vectorstring;
        while (ss >> token)
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

        ss.clear();
    }

    ss.clear();
    // read in the vertices as well as their normals     
    std::string datasentence;
    for (int i=0;i<numvertex;i++)
    {
        std::getline(ifs, datasentence);
        ss.str(datasentence);
        std::string data;
        std::vector<std::string> vectordata;

        vertex v;

        while (ss >> data)
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

        vertices.push_back(v);

        ss.clear();
    }

    ss.clear();
    // read in the triangles
    std::string trianglesentence_;
    for (int i=0;i<numfaces;i++)
    {
        std::getline(ifs,trianglesentence_);
        ss.str(trianglesentence_);
        std::string data;
        std::vector<std::string> vectordata;

        while (ss >> data)
        {
            vectordata.push_back(data);
        }

        ASSERT((vectordata.size() == 4), "The triangle file must contain 3 f1 f2 f3");

        triangle t;

        INT3 faceid;

        for (int j=0;j<3;j++)
        {
            faceid[j] = StringTools::StringToType<int>(vectordata[j+1]);
        }

        t.triangleindices_ = faceid;
        for (int j=0;j<3;j++)
        {
            t.vertices_[j] = vertices[faceid[j]];
        }

        triangles.push_back(t);
        ss.clear();
    }

    ifs.close();

    return true;
}

MeshTools::Real3 MeshTools::calculateShift(const Real3& vec1, const Real3& vec2, const Real3& boxLength)
{
    Real3 diff;
    for (int i=0;i<3;i++)
    {
        Real d = vec1[i] - vec2[i];

        if (d > (0.5 * boxLength[i])) diff[i] = -boxLength[i];
        else if (d < (- 0.5 * boxLength[i])) diff[i] = boxLength[i];
        else diff[i] = 0.0;
    }
    
    return diff;
}

bool MeshTools::isPeriodicEdge(const Real3& vec1, const Real3& vec2, Real3& newarr, const Real3& boxLength)
{
    Real3 diff = calculateShift(vec1, vec2, boxLength);

    bool isPeriodic=false;

    newarr = {};
    for (int i=0;i<3;i++)
    {
        newarr[i] = vec1[i] + diff[i];
    }

    Real3 vecdiff = {};
    for (int i=0;i<3;i++)
    {
        vecdiff[i] = vec1[i] - vec2[i];

        if (std::abs(vecdiff[i]) > (boxLength[i] * 0.5))
        {
            return true;
        }
    }

    return false;
}

void MeshTools::MapVerticesToFaces(Mesh& mesh, std::vector<std::vector<int>>& map)
{
    auto& vertices = mesh.getvertices();
    auto& triangles= mesh.gettriangles();

    map.clear();
    map.resize(vertices.size());

    for (int i=0;i<triangles.size();i++)
    {
        auto t = triangles[i].triangleindices_;

        for (int j=0;j<3;j++)
        {
            map[t[j]].push_back(i);
        }
    }
}

void MeshTools::CalculateTriangleAreasAndFaceNormals(Mesh& mesh, std::vector<Real>& Areas, std::vector<Real3>& Normals)
{
    auto& triangles = mesh.accesstriangles();
    auto& vertices  = mesh.accessvertices();

    Areas.clear();
    Areas.resize(triangles.size());

    Normals.clear();
    Normals.resize(triangles.size());

    // calculate the area of the triangles as well as the normals of the faces of triangles
    for (int i=0;i<triangles.size();i++)
    {
        auto& t = triangles[i];

        int index1 = t.triangleindices_[0];
        int index2 = t.triangleindices_[1];
        int index3 = t.triangleindices_[2];

        Real3 diff1 = {};
        Real3 diff2 = {};
        Real3 diff3 = {};
        Real norm1, norm2, norm3;

        mesh.getVertexDistance(vertices[index1], vertices[index2], diff1, norm1);
        mesh.getVertexDistance(vertices[index3], vertices[index2], diff2, norm2);
        mesh.getVertexDistance(vertices[index3], vertices[index1], diff3, norm3);

        Real3 crossProduct = LinAlg3x3::CrossProduct(diff1, diff2);
        Real norm = LinAlg3x3::norm(crossProduct);
        Real a    = norm*0.5;

        Real3 n;
        for (int j=0;j<3;j++)
        {
            n[j] = crossProduct[j]/norm;
        }

        Normals[i] = n;
        Areas[i]   = a;
    }
}

void MeshTools::CalculateVertexNeighbors(Mesh& mesh, std::vector<std::vector<int>>& neighborIndices)
{
    auto& vertices = mesh.getvertices();
    auto& triangles= mesh.gettriangles();

    neighborIndices.clear();
    neighborIndices.resize(vertices.size());

    for (int i=0;i<triangles.size();i++)
    {
        auto& t = triangles[i];

        for (int j=0;j<3;j++)
        {
            int index1 = t.triangleindices_[j];
            for (int k=0;k<3;k++)
            {
                if (j != k)
                {
                    int index2 = t.triangleindices_[k];

                    auto it = std::find(neighborIndices[index1].begin(), neighborIndices[index1].end(), index2);

                    if (it == neighborIndices[index1].end())
                    {
                        neighborIndices[index1].push_back(index2);
                    }
                }
            }
        }
    }
}

void MeshTools::MapEdgeToFace(Mesh& mesh, std::map<INT2, std::vector<int>>& mapEdgeToFace, \
std::vector<std::vector<INT2>>& mapVertexToEdge)
{
    auto& triangles = mesh.gettriangles();
    auto& vertices  = mesh.getvertices();

    mapVertexToEdge.clear();
    mapVertexToEdge.resize(vertices.size());
    mapEdgeToFace.clear();

    // iterate over all the triangles
    for (int i=0;i<triangles.size();i++)
    {
        // find the triangles indices
        auto& t = triangles[i];

        for (int j=0;j<3;j++)
        {
            int idx1 = t.triangleindices_[j];
            int idx2 = t.triangleindices_[(j+1)%3];
            int minIndex = std::min(idx1, idx2);
            int maxIndex = std::max(idx1, idx2);
            INT2 arr = {{minIndex, maxIndex}};

            auto it = mapEdgeToFace.find(arr);

            if (it == mapEdgeToFace.end())
            {
                std::vector<int> faceIndices;
                faceIndices.push_back(i);

                mapEdgeToFace.insert(std::make_pair(arr,faceIndices));
            }
            else
            {
                it -> second.push_back(i);
            }
        }

        // map vertex to edges 
        for (int j=0;j<3;j++)
        {
            int idx1 = t.triangleindices_[j];
            int idx2 = t.triangleindices_[(j+1)%3];
            int minIdx = std::min(idx1, idx2);
            int maxIdx = std::max(idx1, idx2);

            INT2 indices = {{minIdx, maxIdx}};

            for (int k=0;k<2;k++)
            {
                // find if the edge is already in the vector
                auto f = std::find(mapVertexToEdge[indices[k]].begin(), mapVertexToEdge[indices[k]].end(), indices);

                if (f == mapVertexToEdge[indices[k]].end())
                {
                    mapVertexToEdge[indices[k]].push_back(indices);
                }
            }
        }
    }

    // do a check to make sure that each edge is shared by at least 1 and at most 2 faces
    int id = 0;
    for (auto it=mapEdgeToFace.begin(); it != mapEdgeToFace.end();it++)
    {
        if (it -> second.size() > 2)
        {
            std::cout << "This is for edge " << id << " consisting of vertex " << it ->first[0] <<" " << it->first[1] <<"\n";
            for (auto s : it -> second)
            {
                std::cout << "Face it corresponds to : " << s <<"\n";
            }
        }
        ASSERT((it->second.size() == 1 || it -> second.size() ==2), "The number of faces shared by an edge needs to be either 1 or 2 \
        , however the calculated result shows " << it -> second.size());
        id ++;
    }

    return;
}

void MeshTools::CalculateBoundaryVertices(Mesh& mesh, std::map<INT2, std::vector<int>>& mapEdgeToFace, std::vector<bool>& boundaryIndicator)
{
    const auto& vertices = mesh.getvertices();

    boundaryIndicator.clear();
    boundaryIndicator.resize(vertices.size(), false);

    for (auto it = mapEdgeToFace.begin(); it != mapEdgeToFace.end(); it ++)
    {
        int size = it -> second.size();
        ASSERT(( size == 1 || size == 2), "An edge can only be shared by 1 or 2 faces.");
        INT2 arr = it -> first;

        // if the edge size == 1, then it is a boundary edge
        if (size == 1)
        {
            boundaryIndicator[arr[0]] = true;
            boundaryIndicator[arr[1]] = true;
        }
    }
}

bool MeshTools::IsBoundary(int index, const std::vector<bool>& boundaryIndicator)
{
    return boundaryIndicator[index];
}

void MeshTools::ConvertToNonPBCMesh(Mesh& mesh, std::vector<Real3>& vertices, std::vector<INT3>& triangles)
{
    // if periodic, then do something, else do nothing 
    if (mesh.isPeriodic())
    {
        auto MeshVertices = mesh.getvertices();
        auto MeshTriangles = mesh.gettriangles();

        // first let's copy all the vertices 
        vertices.clear();
        for (auto& v : MeshVertices)
        {
            vertices.push_back(v.position_);
        }

        // check if triangle is periodic 
        for (auto& t : MeshTriangles)
        {
            // find all the edge lengths
            bool periodicTriangle = MeshTools::IsPeriodicTriangle(MeshVertices, t.triangleindices_, mesh.getBoxLength());

            // if this particular triangle is not periodic 
            if (! periodicTriangle)
            {
                triangles.push_back(t.triangleindices_);
            }
            // if it's periodic triangle, then we push back 3 new vertices 
            else
            {
                Real3 verticesNew1;
                Real3 verticesNew2, verticesDiff2;
                Real distsq2;
                Real3 verticesNew3, verticesDiff3;
                Real distsq3;
                int idx1 = t.triangleindices_[0];
                int idx2 = t.triangleindices_[1];
                int idx3 = t.triangleindices_[2];

                // get the pbc corrected distance 
                mesh.getVertexDistance(MeshVertices[idx2].position_, MeshVertices[idx1].position_,verticesDiff2, distsq2);
                mesh.getVertexDistance(MeshVertices[idx3].position_, MeshVertices[idx1].position_,verticesDiff3, distsq3); 

                // get the new vertices --> with respect to position 1
                for (int j=0;j<3;j++)
                {
                    verticesNew2[j] = MeshVertices[idx1].position_[j] + verticesDiff2[j];
                    verticesNew3[j] = MeshVertices[idx1].position_[j] + verticesDiff3[j];
                }

                // Find approximately the center of the triangle
                Real3 center_of_triangle = {};
                for (int j=0;j<3;j++)
                {
                    center_of_triangle[j] += MeshVertices[idx1].position_[j];
                    center_of_triangle[j] += verticesNew2[j];
                    center_of_triangle[j] += verticesNew3[j];
                }

                for (int j=0;j<3;j++)
                {
                    center_of_triangle[j] *= 1.0/3.0;
                }

                Real3 shift = mesh.getShiftIntoBox(center_of_triangle);

                for (int j=0;j<3;j++)
                {
                    verticesNew1[j] = MeshVertices[idx1].position_[j] + shift[j];
                    verticesNew2[j] = verticesNew2[j] + shift[j];
                    verticesNew3[j] = verticesNew3[j] + shift[j];
                }

                int NewIndex1 = vertices.size();
                vertices.push_back(verticesNew1);
                int NewIndex2 = vertices.size();
                vertices.push_back(verticesNew2);
                int NewIndex3 = vertices.size();
                vertices.push_back(verticesNew3);

                INT3 NewT = {{NewIndex1, NewIndex2, NewIndex3}};
                triangles.push_back(NewT);
            }
        }
    }
}

bool MeshTools::IsPeriodicTriangle(std::vector<vertex>& Vertices,INT3& face, Real3 BoxLength)
{
    bool IsPeriodic=false;
    for (int i=0;i<3;i++)
    {
        int index1 = face[i];
        int index2 = face[(i+1) % 3]; 
        Real3 diff;
        for (int j=0;j<3;j++)
        {
            diff[j] = Vertices[index1].position_[j] - Vertices[index2].position_[j];

            if (std::abs(diff[j]) >= 0.5 * BoxLength[j])
            {
                return true;
            }
        }
    }

    return false;
}

MeshTools::INT2 MeshTools::makeEdge(int i, int j)
{
    int minIndex = std::min(i,j);
    int maxIndex = std::max(i,j);

    INT2 ret = {{minIndex, maxIndex}};

    return ret;
}

void MeshTools::MapEdgeToOpposingVertices(Mesh& mesh, std::map<INT2, std::vector<int>>& mapEdgeToFace, std::map<INT2, std::vector<int>>& MapEdgeToOppoVertices)
{
    const auto& tri = mesh.gettriangles();

    // clear the map at which we want to write to  
    MapEdgeToOppoVertices.clear();

    // iterate over the edges 
    for (auto it = mapEdgeToFace.begin(); it != mapEdgeToFace.end(); it ++)
    {
        std::vector<int> faces = it -> second;
        INT2 edge = it -> first;
        std::vector<int> opposingPoints;

        int numf = faces.size();

        // only perform calculation if we are working with non boundary edge --> (edges shared by 2 faces)
        if (numf == 2)
        {
            // iterate over the 2 faces 
            for (int f : faces)
            {
                auto& TriIndices = tri[f].triangleindices_;

                // iterate over the triangular indices of each of the face 
                for (int id : TriIndices)
                {
                    bool IsInEdge = std::find(edge.begin(), edge.end(), id) != edge.end();

                    if ( ! IsInEdge)
                    {
                        opposingPoints.push_back(id);
                    }
                }
            }

            ASSERT((opposingPoints.size() == 2), "The opposing points of an edge needs to be 2 while it is " << opposingPoints.size());
            MapEdgeToOppoVertices.insert(std::make_pair(edge, opposingPoints));
        }
    }
}

void MeshTools::CalculateCornerArea(Mesh& mesh, std::vector<Real3>& CornerArea, std::vector<Real>& VertexArea)
{
    const auto& vertices = mesh.getvertices();
    const auto& triangles= mesh.gettriangles();
    int nf = triangles.size();
    int nv = vertices.size();

    // resize corner area
    CornerArea.clear();
    CornerArea.resize(nf);

    // resize vertex area
    VertexArea.clear();
    VertexArea.resize(nv,0.0);

	for (int i = 0; i < nf; i++) {
		// Edges
        auto& t = triangles[i];
        int index1 = t.triangleindices_[0];
        int index2 = t.triangleindices_[1];
        int index3 = t.triangleindices_[2];

        Real3 edge1, edge2, edge3;
        Real edge1sq, edge2sq, edge3sq;

        mesh.getVertexDistance(vertices[index3], vertices[index2], edge1, edge1sq);
        mesh.getVertexDistance(vertices[index1], vertices[index3], edge2, edge2sq);
        mesh.getVertexDistance(vertices[index2], vertices[index1], edge3, edge3sq);

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
			CornerArea[i][1] = -0.25f * l2[2] * area/LinAlg3x3::DotProduct(edge1, edge3);
			CornerArea[i][2] = -0.25f * l2[1] * area/LinAlg3x3::DotProduct(edge1, edge2);
			CornerArea[i][0] = area - CornerArea[i][1] - CornerArea[i][2];
		} else if (bcw[1] <= 0.0f) {
			CornerArea[i][2] = -0.25f * l2[0] * area/LinAlg3x3::DotProduct(edge2, edge1);
			CornerArea[i][0] = -0.25f * l2[2] * area/LinAlg3x3::DotProduct(edge2, edge3);
			CornerArea[i][1] = area - CornerArea[i][2] - CornerArea[i][0];
		} else if (bcw[2] <= 0.0f) {
			CornerArea[i][0] = -0.25f * l2[1] * area/LinAlg3x3::DotProduct(edge3, edge2);
			CornerArea[i][1] = -0.25f * l2[0] * area/LinAlg3x3::DotProduct(edge3, edge1);
			CornerArea[i][2] = area - CornerArea[i][0] - CornerArea[i][1];
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

				CornerArea[i][j] = scale * (bcw[next] +
				                             bcw[nextnext]);
            }
		}

		VertexArea[t.triangleindices_[0]] += CornerArea[i][0];
		VertexArea[t.triangleindices_[1]] += CornerArea[i][1];
		VertexArea[t.triangleindices_[2]] += CornerArea[i][2];
	}
}
