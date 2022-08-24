#include "MeshActions.h"

void MeshActions::TranslateMesh(CommandLineArguments& cmd)
{
    std::string inputfname;
    std::string outputfname;

    Mesh mesh;
    std::vector<Real> translate_vec;
    cmd.readString("i", CommandLineArguments::Keys::Required, inputfname);
    cmd.readString("o", CommandLineArguments::Keys::Optional, outputfname);
    cmd.readVector("trans", CommandLineArguments::Keys::Required, translate_vec);

    ASSERT((translate_vec.size()==3), "The inputted translation vector must be of size 3");

    MeshTools::readPLYlibr(inputfname, mesh);
    std::vector<Real3> vertices;
    std::vector<INT3> faces;
    std::vector<Real3> normals;

    const auto& verts = mesh.getvertices();
    const auto& tri   = mesh.gettriangles();

    vertices.resize(verts.size());
    faces.resize(tri.size());
    normals.resize(verts.size());

    #pragma omp parallel for
    for (int i=0;i<vertices.size();i++)
    {
        for (int j=0;j<3;j++)
        {
            vertices[i][j] = verts[i].position_[j] + translate_vec[j];
            normals[i][j]  = verts[i].normals_[j];
        }
    }

    for(int i=0;i<faces.size();i++)
    {
        faces[i] = tri[i].triangleindices_;
    }

    MeshTools::writePLY(outputfname, vertices, faces, normals);
}

void MeshActions::CurveFit(CommandLineArguments& cmd)
{
    std::string inputfname;
    std::string outputfname="curvefit.out";
    std::string faceCfname="curvefit_faceC.out";
    std::string ff2fname="ff2.out";
    std::string neighbors;
    ParameterPack pack;

    Mesh mesh;
    Real3 box;
    curveptr curve;

    cmd.readString("i", CommandLineArguments::Keys::Required, inputfname);
    cmd.readString("o", CommandLineArguments::Keys::Optional, outputfname);
    bool fc_read  = cmd.readString("fc", CommandLineArguments::Keys::Optional, faceCfname);
    bool ff2_read = cmd.readString("ff2", CommandLineArguments::Keys::Optional, ff2fname);
    cmd.readString("neighbors", CommandLineArguments::Keys::Required, neighbors);
    bool pbcMesh = cmd.readArray("box", CommandLineArguments::Keys::Optional, box);

    // insert the values into parameter pack
    pack.insert("neighbors", neighbors);

    if (pbcMesh)
    {
        mesh.setBoxLength(box);
    }

    // read the mesh
    MeshTools::readPLYlibr(inputfname, mesh);

    // initialize curvature
    CurvatureInput input = {{pack}};
    curve = curveptr(CurvatureRegistry::Factory::instance().create("curvefit", input));

    // calculate the curvature 
    curve->calculate(mesh);

    // print the curvature 
    curve->printCurvature(outputfname);

    if (fc_read)
    {
        curve->printFaceCurvature(faceCfname);
    }

    if (ff2_read)
    {
        CurvatureCurveFit* c = dynamic_cast<CurvatureCurveFit*>(curve.get());
        ASSERT((c != nullptr), "Something went wrong in curvefit mesh action.");
        c->printff2(ff2fname);
    }

}

void MeshActions::JetFit(CommandLineArguments& cmd)
{
    std::string inputfname;
    std::string outputfname="jetfit.out";
    std::string faceCfname="jetfit_faceC.out";
    std::string neighbors;
    std::string degree;
    std::string MongeCoefficient;

    ParameterPack pack;

    Mesh mesh;
    Real3 box;
    curveptr curve;

    cmd.readString("i", CommandLineArguments::Keys::Required, inputfname);
    cmd.readString("o", CommandLineArguments::Keys::Optional, outputfname);
    cmd.readString("fc", CommandLineArguments::Keys::Optional, faceCfname);
    cmd.readString("neighbors", CommandLineArguments::Keys::Required, neighbors);
    cmd.readString("degree", CommandLineArguments::Keys::Required, degree);
    cmd.readString("MongeCoefficient", CommandLineArguments::Keys::Required, MongeCoefficient);
    bool pbcMesh = cmd.readArray("box", CommandLineArguments::Keys::Optional, box);

    // insert the values into parameter pack
    pack.insert("neighbors", neighbors);
    pack.insert("degree", degree);
    pack.insert("MongeCoefficient", MongeCoefficient);

    if (pbcMesh)
    {
        mesh.setBoxLength(box);
    }

    // read the mesh
    MeshTools::readPLYlibr(inputfname, mesh);

    // initialize curvature
    CurvatureInput input = {{pack}};
    curve = curveptr(CurvatureRegistry::Factory::instance().create("jetfit", input));

    curve->calculate(mesh);

    curve->printCurvature(outputfname);
    curve->printFaceCurvature(faceCfname);
}

void MeshActions::CurvatureFlow(CommandLineArguments& cmd)
{
    refineptr refine;
    std::string inputfname;
    std::string outputfname="CurvatureFlow.ply";
    std::string iterations;
    std::string lambdadt="1";
    Real3 box;
    ParameterPack pack;
    bool nonpbc = true;
    Mesh mesh;

    cmd.readString("i", CommandLineArguments::Keys::Required, inputfname);
    cmd.readString("o", CommandLineArguments::Keys::Optional, outputfname);
    cmd.readString("iteration", CommandLineArguments::Keys::Required, iterations);
    cmd.readString("lambdadt", CommandLineArguments::Keys::Optional, lambdadt);
    cmd.readBool("nonpbc", CommandLineArguments::Keys::Optional, nonpbc);
    bool pbcMesh = cmd.readArray("box", CommandLineArguments::Keys::Optional, box);
    if (pbcMesh)
    {
        mesh.setBoxLength(box);
    }

    // read the mesh
    MeshTools::readPLYlibr(inputfname, mesh);

    // insert the values into parameter pack
    pack.insert("iterations", iterations);
    pack.insert("lambdadt", lambdadt);
    pack.insert("name", "temp");

    MeshRefineStrategyInput input = {{pack}};
    refine = refineptr(MeshRefineStrategyFactory::Factory::instance().create("curvatureflow", input));

    refine -> refine(mesh);

    if (nonpbc)
    {
        mesh.printNonPBCMesh(outputfname);
    }
    else
    {
        mesh.printPLY(outputfname);
    }
}

void MeshActions::FindNonPBCTriangles(CommandLineArguments& cmd)
{
    std::string inputfname;
    std::string outputfname="nonpbc.out";
    Real3 box;
    ParameterPack pack;
    Mesh mesh;

    cmd.readString("i", CommandLineArguments::Keys::Required, inputfname);
    cmd.readString("o", CommandLineArguments::Keys::Optional, outputfname);
    bool pbcMesh = cmd.readArray("box", CommandLineArguments::Keys::Required, box);
    if (pbcMesh)
    {
        mesh.setBoxLength(box);
    }

    // read the mesh
    MeshTools::readPLYlibr(inputfname, mesh);

    mesh.printNonPeriodicTriangleIndices(outputfname);
}

void MeshActions::CutMesh(CommandLineArguments& cmd)
{
    std::string inputfname;
    std::string outputfname="out.ply";
    Mesh mesh;
    Real3 volume;

    cmd.readString("i", CommandLineArguments::Keys::Required, inputfname);
    cmd.readString("o", CommandLineArguments::Keys::Optional, outputfname);
    cmd.readArray("volume", CommandLineArguments::Keys::Required, volume);


    // read the mesh
    MeshTools::readPLYlibr(inputfname, mesh);

    std::vector<INT3> face;
    std::vector<Real3> verts;
    MeshTools::CutMesh(mesh, face, verts, volume);
    MeshTools::writePLY(outputfname, verts, face);
}

void MeshActions::ConvertToNonPBCMesh(CommandLineArguments& cmd)
{
    std::string inputfname;
    std::string outputfname="nonpbc.ply";

    Mesh mesh;
    Real3 Box;

    cmd.readString("i", CommandLineArguments::Keys::Required, inputfname);
    cmd.readString("o", CommandLineArguments::Keys::Optional, outputfname);
    cmd.readArray("box", CommandLineArguments::Keys::Required, Box);

    mesh.setBoxLength(Box);

    // read the mesh 
    MeshTools::readPLYlibr(inputfname, mesh);

    // output non pbc mesh 
    std::vector<INT3> face;
    std::vector<Real3> verts;
    MeshTools::ConvertToNonPBCMesh(mesh, verts, face);

    // write the non pbc mesh 
    MeshTools::writePLY(outputfname, verts, face);
}

void MeshActions::ScaleMesh(CommandLineArguments& cmd)
{
    std::string inputfname;
    std::string outputfname="scaled.ply";

    Mesh mesh;
    Real scale;

    cmd.readString("i", CommandLineArguments::Keys::Required, inputfname);
    cmd.readString("o", CommandLineArguments::Keys::Optional, outputfname);
    cmd.readValue("scale", CommandLineArguments::Keys::Required, scale);

    std::vector<INT3> face;
    std::vector<Real3> verts;

    // read the Mesh
    MeshTools::readPLYlibr(inputfname, mesh);

    const auto& meshVerts = mesh.getvertices();
    const auto& meshFaces = mesh.gettriangles();

    for (int i=0;i<meshVerts.size();i++)
    {
        Real3 pos;
        for (int j=0;j<3;j++)
        {
            pos[j] = meshVerts[i].position_[j] * scale;
        }
        verts.push_back(pos);
    }

    for (int i=0;i<meshFaces.size();i++)
    {
        face.push_back(meshFaces[i].triangleindices_);
    }

    // write to output
    MeshTools::writePLY(outputfname, verts, face);
}

void MeshActions::ColorVertex(CommandLineArguments& cmd)
{

}

void MeshActions::Project3dCurvature(CommandLineArguments& cmd)
{
    using INT2 = std::array<int,2>;

    std::string inputfname, outputfname, fcfname;
    Real origin_pos=0.0;
    int index=0;
    int n1, n2, colnum;
    Real L1, L2, d1, d2;
    std::vector<Real> fc;
    Real3 box;

    cmd.readValue("i", CommandLineArguments::Keys::Required, inputfname);
    cmd.readValue("o", CommandLineArguments::Keys::Optional, outputfname);
    cmd.readValue("index", CommandLineArguments::Keys::Optional, index);
    cmd.readValue("origin", CommandLineArguments::Keys::Optional, origin_pos);
    cmd.readValue("n1", CommandLineArguments::Keys::Required, n1);
    cmd.readValue("n2", CommandLineArguments::Keys::Required, n2);
    cmd.readValue("L1", CommandLineArguments::Keys::Required, L1);
    cmd.readValue("L2", CommandLineArguments::Keys::Required, L2);
    cmd.readValue("fc", CommandLineArguments::Keys::Required, fcfname);
    cmd.readValue("col", CommandLineArguments::Keys::Required, colnum);
    bool isPBC = cmd.readArray("box", CommandLineArguments::Keys::Optional, box);

    // read the tabulated data for curvature
    StringTools::ReadTabulatedData(fcfname, colnum,fc);

    // read the input file and set up the mesh 
    Mesh mesh;
    MeshTools::readPLYlibr(inputfname, mesh);
    if (isPBC)
    {
        mesh.setBoxLength(box);
    }

    // direction 
    Real3 D = {{0,0,0}};
    D[index] = 1;

    // other index 
    std::vector<int> otherIndex;
    for (int i=0;i<3;i++)
    {
        if (i != index)
        {
            otherIndex.push_back(i);
        }
    }

    // calculate the delta 
    d1 = L1 / n1;
    d2 = L2 / n2;

    // initialize the points
    std::vector<Real3> points(n1*n2);
    std::vector<Real> pointCurvature(n1*n2);

    // populate the points 
    for (int i=0;i<n1;i++)
    {
        for (int j=0;j<n2;j++)
        {
            Real3 p;
            p[index] = origin_pos;
            p[otherIndex[0]] = i * d1;
            p[otherIndex[1]] = j * d2;

            points.push_back(p);
        }
    }

    // faces, verts of the mesh 
    const auto& face = mesh.gettriangles();
    const auto& verts= mesh.getvertices();
    const auto& faceNormal = mesh.getFaceNormals();

    #pragma omp parallel for 
    for (int i=0;i<points.size();i++)
    {
        Real3 O = points[i];
        std::vector<int> faceIndex;
        std::vector<Real> value;
        for (int j=0;j<face.size();j++)
        {
            Real3 A = verts[face[j].triangleindices_[0]].position_;
            Real3 B, C;
            Real t, u, v;
            if (isPBC)
            {
                B = mesh.getShiftedVertexPosition(verts[face[j].triangleindices_[0]], verts[face[j].triangleindices_[1]]);
                C = mesh.getShiftedVertexPosition(verts[face[j].triangleindices_[0]], verts[face[j].triangleindices_[2]]);
            }
            else
            {
                B = verts[face[j].triangleindices_[1]].position_;
                C = verts[face[j].triangleindices_[2]].position_;
            }

            bool isIntersect = MeshTools::MTRayTriangleIntersection(A, B, C, \
                                                O, D, faceNormal[j], \
                                                t, u ,v);
            if (isIntersect)
            {
                faceIndex.push_back(j);
                value.push_back(t);
            }
        }

        // find the min element of value
        auto it = std::min_element(value.begin(), value.end());
        int minIndex = it - value.begin();
        pointCurvature[i] = fc[minIndex];
    }
}