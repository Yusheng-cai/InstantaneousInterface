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
    std::string neighbors;
    ParameterPack pack;

    Mesh mesh;
    Real3 box;
    curveptr curve;

    cmd.readString("i", CommandLineArguments::Keys::Required, inputfname);
    cmd.readString("o", CommandLineArguments::Keys::Optional, outputfname);
    cmd.readString("fc", CommandLineArguments::Keys::Optional, faceCfname);
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

    curve->calculate(mesh);

    curve->printCurvature(outputfname);
    curve->printFaceCurvature(faceCfname);
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
    cmd.readString("iterations", CommandLineArguments::Keys::Required, iterations);
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