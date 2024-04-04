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
    for (int i=0;i<vertices.size();i++){
        for (int j=0;j<3;j++){
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

void MeshActions::QuadraticCurveFit(CommandLineArguments& cmd){
    std::string inputfname;
    std::string outputfname="quadraticfit.out";
    std::string neighbors, MonteCarlo="true", MonteCarloN="50";
    ParameterPack pack;

    Mesh mesh;
    Real3 box;
    curveptr curve;

    cmd.readString("i", CommandLineArguments::Keys::Required, inputfname);
    cmd.readString("o", CommandLineArguments::Keys::Optional, outputfname);
    cmd.readString("neighbors", CommandLineArguments::Keys::Required, neighbors);
    cmd.readString("MonteCarlo", CommandLineArguments::Keys::Optional, MonteCarlo);
    cmd.readString("MonteCarloN", CommandLineArguments::Keys::Optional, MonteCarloN);
    bool pbcMesh = cmd.readArray("box", CommandLineArguments::Keys::Optional, box);

    // insert the values into parameter pack
    pack.insert("neighbors", neighbors);
    pack.insert("MonteCarlo", MonteCarlo);
    pack.insert("MonteCarloN", MonteCarloN);

    // read the mesh
    MeshTools::readPLYlibr(inputfname, mesh);
    if (pbcMesh){
        mesh.setBoxLength(box);
    }

    // initialize curvature
    CurvatureInput input = {{pack}};
    curve = curveptr(CurvatureRegistry::Factory::instance().create("quadraticfit", input));

    // calculate the curvature 
    curve->calculate(mesh);

    // print the curvature 
    curve->printCurvature(outputfname);
}

void MeshActions::CurveFit(CommandLineArguments& cmd)
{
    std::string inputfname;
    std::string outputfname="curvefit.out";
    std::string faceCfname="curvefit_faceC.out";
    std::string ff2fname="ff2.out";
    std::string neighbors, boundaryNeighbors;
    ParameterPack pack;

    Mesh mesh;
    Real3 box;
    curveptr curve;
    bool extend_boundary=false;

    cmd.readString("i", CommandLineArguments::Keys::Required, inputfname);
    cmd.readString("o", CommandLineArguments::Keys::Optional, outputfname);
    extend_boundary = cmd.readString("boundaryNeighbor", CommandLineArguments::Keys::Optional, boundaryNeighbors);
    bool fc_read  = cmd.readString("fc", CommandLineArguments::Keys::Optional, faceCfname);
    bool ff2_read = cmd.readString("ff2", CommandLineArguments::Keys::Optional, ff2fname);
    cmd.readString("neighbors", CommandLineArguments::Keys::Required, neighbors);
    bool pbcMesh = cmd.readArray("box", CommandLineArguments::Keys::Optional, box);

    // insert the values into parameter pack
    pack.insert("neighbors", neighbors);
    if (extend_boundary){
        pack.insert("extend_boundary", "true");
        pack.insert("boundary_neighbors", boundaryNeighbors);
    }

    if (pbcMesh){
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

    if (fc_read){
        curve->printFaceCurvature(faceCfname);
    }

    if (ff2_read){
        CurvatureCurveFit* c = dynamic_cast<CurvatureCurveFit*>(curve.get());
        ASSERT((c != nullptr), "Something went wrong in curvefit mesh action.");
        c-> printSecondFundamentalForm(ff2fname);
    }
}

void MeshActions::TensorFit(CommandLineArguments& cmd)
{
    std::string inputfname, outputfname="tensor.out", faceCfname="tensor_faceC.out";

    Real3 box;

    // Curvature parameters pack
    ParameterPack pack;

    // initialize the necessary pointers
    Mesh mesh;
    curveptr curve;

    cmd.readValue("i", CommandLineArguments::Keys::Required, inputfname);
    cmd.readValue("o", CommandLineArguments::Keys::Optional, outputfname);
    bool pbcMesh = cmd.readArray("box", CommandLineArguments::Keys::Optional, box);
    bool fcread  = cmd.readString("fc", CommandLineArguments::Keys::Optional, faceCfname);
    if (pbcMesh){
        mesh.setBoxLength(box);
    }

    // read the mesh
    MeshTools::readPLYlibr(inputfname, mesh);

    // initialize curvature
    CurvatureInput input = {{pack}};
    curve = curveptr(CurvatureRegistry::Factory::instance().create("tensor", input));

    curve->calculate(mesh);

    curve->printCurvature(outputfname);

    if (fcread){
        curve->printFaceCurvature(faceCfname);
    }
}

void MeshActions::JetFit(CommandLineArguments& cmd)
{
    std::string inputfname, outputfname="jetfit.out", faceCfname="jetfit_faceC.out", neighbors, degree, MongeCoefficient;
    std::string pdir_fname="pdir.out", pca_fname;
    std::string monge;
    Real3 box;

    // Curvature parameters pack
    ParameterPack pack;

    // initialize the necessary pointers
    Mesh mesh;
    curveptr curve;

    cmd.readString("i", CommandLineArguments::Keys::Required, inputfname);
    cmd.readString("o", CommandLineArguments::Keys::Optional, outputfname);
    bool fcread = cmd.readString("fc", CommandLineArguments::Keys::Optional, faceCfname);
    bool MongeRead = cmd.readString("mongeFile", CommandLineArguments::Keys::Optional, monge);
    bool pdirRead  = cmd.readString("pdirFile", CommandLineArguments::Keys::Optional, pdir_fname);
    bool pcaRead   = cmd.readString("PCAFile", CommandLineArguments::Keys::Optional, pca_fname);
    cmd.readString("neighbors", CommandLineArguments::Keys::Required, neighbors);
    cmd.readString("degree", CommandLineArguments::Keys::Required, degree);
    cmd.readString("MongeCoefficient", CommandLineArguments::Keys::Required, MongeCoefficient);
    bool pbcMesh = cmd.readArray("box", CommandLineArguments::Keys::Optional, box);

    // insert the values into parameter pack
    pack.insert("neighbors", neighbors);
    pack.insert("degree", degree);
    pack.insert("MongeCoefficient", MongeCoefficient);

    if (pbcMesh){
        mesh.setBoxLength(box);
    }

    // read the mesh
    MeshTools::readPLYlibr(inputfname, mesh);

    // initialize curvature
    CurvatureInput input = {{pack}};
    curve = curveptr(CurvatureRegistry::Factory::instance().create("jetfit", input));

    curve->calculate(mesh);

    curve->printCurvature(outputfname);

    if (fcread){
        curve->printFaceCurvature(faceCfname);
    }

    if (MongeRead){
        CurvatureJetFit* jetfit = dynamic_cast<CurvatureJetFit*>(curve.get());
        ASSERT((jetfit != nullptr), "Something went wrong in jetfit.");
        jetfit->printCoefficientPerVertex(monge);
    }

    if (pcaRead){
        CurvatureJetFit* jetfit = dynamic_cast<CurvatureJetFit*>(curve.get());
        ASSERT((jetfit != nullptr), "Something went wrong in jetfit.");
        jetfit->printPCAeigenvector(pca_fname);
    }

    if (pdirRead){
        curve -> printPrincipalDir(pdir_fname);
    }
}

void MeshActions::FDMFit(CommandLineArguments& cmd){
    std::string inputfname, outputfname="FDMfit.out", faceCfname="FDMfit_faceC.out";
    Real3 box;

    // Curvature parameters pack
    ParameterPack pack;

    // initialize the necessary pointers
    Mesh mesh;
    curveptr curve;

    cmd.readString("i", CommandLineArguments::Keys::Required, inputfname);
    cmd.readString("o", CommandLineArguments::Keys::Optional, outputfname);
    bool fcread = cmd.readString("fc", CommandLineArguments::Keys::Optional, faceCfname);
    bool pbcMesh = cmd.readArray("box", CommandLineArguments::Keys::Optional, box);

    if (pbcMesh){
        mesh.setBoxLength(box);
    }

    // read the mesh
    MeshTools::readPLYlibr(inputfname, mesh);

    // initialize curvature
    CurvatureInput input = {{pack}};
    curve = curveptr(CurvatureRegistry::Factory::instance().create("FDM", input));

    curve->calculate(mesh);

    curve->printCurvature(outputfname);

    if (fcread)
    {
        curve->printFaceCurvature(faceCfname);
    }
}

void MeshActions::CurvatureFlow(CommandLineArguments& cmd)
{
    refineptr refine;
    std::string inputfname, outputfname="CurvatureFlow.ply", fixed_index_file, iterations, lambdadt, decimate="true", numBoundarySmooth="0", scale="false", shift_COM="false";
    Real3 box;
    ParameterPack pack;
    bool pbcOutput = true;

    cmd.readString("i", CommandLineArguments::Keys::Required, inputfname);
    cmd.readString("o", CommandLineArguments::Keys::Optional, outputfname);
    cmd.readString("iteration", CommandLineArguments::Keys::Required, iterations);
    cmd.readString("lambdadt", CommandLineArguments::Keys::Required, lambdadt);
    cmd.readBool("pbcOutput", CommandLineArguments::Keys::Optional, pbcOutput);
    cmd.readString("Decimate", CommandLineArguments::Keys::Optional, decimate);
    cmd.readString("NumBoundarySmoothing", CommandLineArguments::Keys::Optional, numBoundarySmooth);
    cmd.readString("scale", CommandLineArguments::Keys::Optional, scale);
    cmd.readString("shiftCOM", CommandLineArguments::Keys::Optional, shift_COM);
    if (cmd.readString("fixed_index_file", CommandLineArguments::Keys::Optional, fixed_index_file)){
        pack.insert("fixed_index_file", fixed_index_file);
    }

    // check if the mesh is periodic 
    bool pbcMesh = cmd.readArray("box", CommandLineArguments::Keys::Optional, box);

    // read the mesh
    Mesh mesh;
    MeshTools::readPLYlibr(inputfname, mesh);
    if (pbcMesh){mesh.setBoxLength(box);}

    // insert the values into parameter pack
    pack.insert("iterations", iterations);
    pack.insert("lambdadt", lambdadt);
    pack.insert("name", "temp");
    pack.insert("Decimate", decimate);
    pack.insert("NumBoundarySmoothing",numBoundarySmooth);
    pack.insert("scale", scale);
    pack.insert("shiftCOM", shift_COM);


    MeshRefineStrategyInput input = {{pack}};
    refine = refineptr(MeshRefineStrategyFactory::Factory::instance().create("curvatureflow", input));

    // refine the mesh
    refine -> refine(mesh);

    if (pbcMesh){
        if (pbcOutput){
            MeshTools::writePLY(outputfname, mesh);
        }
        else {MeshTools::writeNonPBCMesh(outputfname, mesh);}
    }
    else{MeshTools::writePLY(outputfname, mesh);}
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
    if (pbcMesh){
        mesh.setBoxLength(box);
    }

    // read the mesh
    MeshTools::readPLYlibr(inputfname, mesh);

    MeshTools::writeNonPeriodicTriangleIndices(outputfname, mesh);
}

void MeshActions::CutMesh(CommandLineArguments& cmd){
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

    bool addnewtriangles=true;

    cmd.readString("i", CommandLineArguments::Keys::Required, inputfname);
    cmd.readString("o", CommandLineArguments::Keys::Optional, outputfname);
    cmd.readArray("box", CommandLineArguments::Keys::Required, Box);
    cmd.readBool("AddNewTriangle", CommandLineArguments::Keys::Optional, addnewtriangles);


    mesh.setBoxLength(Box);

    // read the mesh 
    MeshTools::readPLYlibr(inputfname, mesh);

    // output non pbc mesh 
    std::vector<INT3> face;
    std::vector<Real3> verts;
    MeshTools::ConvertToNonPBCMesh(mesh, addnewtriangles);
    mesh.CalcVertexNormals();

    // write the non pbc mesh 
    // std::string ext = StringTools::ReadFileExtension(outputfname);
    MeshTools::writePLY(outputfname, mesh);
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

    for (int i=0;i<meshVerts.size();i++){
        Real3 pos;
        for (int j=0;j<3;j++){
            pos[j] = meshVerts[i].position_[j] * scale;
        }
        verts.push_back(pos);
    }

    for (int i=0;i<meshFaces.size();i++){face.push_back(meshFaces[i].triangleindices_);}

    // write to output
    MeshTools::writePLY(outputfname, verts, face);
}

void MeshActions::ColorVertex(CommandLineArguments& cmd)
{
    std::string inputfname, cfname, outputfname;
    int col;
    Real vmin, vmax;
    Real3 box;

    // read in the necessary inputs 
    cmd.readValue("curvature",CommandLineArguments::Keys::Required, cfname);
    cmd.readValue("i", CommandLineArguments::Keys::Required, inputfname);
    cmd.readValue("o", CommandLineArguments::Keys::Optional, outputfname);
    cmd.readValue("col", CommandLineArguments::Keys::Required, col);
    bool isPBC = cmd.readArray("box", CommandLineArguments::Keys::Optional, box);
    bool vminRead = cmd.readValue("vmin", CommandLineArguments::Keys::Optional, vmin);
    bool vmaxRead = cmd.readValue("vmax", CommandLineArguments::Keys::Optional, vmax);

    // read in the curvature 
    std::vector<Real> curvature;
    StringTools::ReadTabulatedData(cfname, col, curvature);

    // read in the mesh 
    Mesh mesh;
    MeshTools::readPLYlibr(inputfname, mesh);

    // find the min and max curvature 
    if (! vminRead)
    {
        vmin = *std::min_element(curvature.begin(), curvature.end());
    }

    if (! vmaxRead)
    {
        vmax = *std::max_element(curvature.begin(), curvature.end());
    }

    // initialize RGB value
    std::vector<Real3> RGB(curvature.size());
    for (int i=0;i<curvature.size();i++)
    {
        Real3 c={{0,0,0}};
        if (curvature[i] > vmax || curvature[i] < vmin)
        {
            c = {{0,0,0}};
        }
        else
        {
            c[0] = 255.0;
            c[1] = 0.0;
            c[2] = std::round(255 * (curvature[i] - vmin)/(vmax - vmin));
        }

        RGB[i] = c;
    }

    // copy vertices and faces 
    std::vector<Real3> v;
    std::vector<INT3> f;

    // nonpbc index
    const auto& verts = mesh.getvertices();
    const auto& faces = mesh.gettriangles();
    // the triangles which are non pbc 
    std::vector<int> nonpbcIndex;
    if (isPBC)
    {
        for (int i=0;i<faces.size();i++)
        {
            if (! MeshTools::IsPeriodicTriangle(const_cast<std::vector<vertex>&>(verts), const_cast<INT3&>(faces[i].triangleindices_),\
                                            box))
            {
                nonpbcIndex.push_back(i);

            }
        }
    }
    else
    {
        nonpbcIndex.resize(faces.size());
        std::iota(nonpbcIndex.begin(), nonpbcIndex.end(), 0);
    }
    for (int i=0;i<verts.size();i++)
    {
        v.push_back(verts[i].position_);
    }
    for (int i=0;i<nonpbcIndex.size();i++)
    {
        f.push_back(faces[nonpbcIndex[i]].triangleindices_);
    }
    
    MeshTools::writePLYRGB(outputfname, v, f, RGB);
}

void MeshActions::Project3dMesh(CommandLineArguments& cmd)
{
    using INT2 = std::array<int,2>;

    // define input output file names 
    std::vector<std::string> inputfname;
    std::string outputfname;
    std::string normalMapfname;

    // some large number 
    Real height_pos=100.0;
    int colnum;
    Real2 L, n, d, ProjectedIndex;
    std::vector<Real> fc;
    Real3 box, D;

    // read the inputs 
    cmd.readVector("i", CommandLineArguments::Keys::Required, inputfname);
    cmd.readValue("o", CommandLineArguments::Keys::Optional, outputfname);
    cmd.readArray("RayDirection", CommandLineArguments::Keys::Required,D);
    LinAlg3x3::normalize(D);
    cmd.readArray("ProjectedIndex", CommandLineArguments::Keys::Required, ProjectedIndex);
    bool normalMap = cmd.readValue("normalMapOutput", CommandLineArguments::Keys::Optional, normalMapfname);

    // the origin 
    cmd.readValue("height", CommandLineArguments::Keys::Optional, height_pos);
    cmd.readArray("L", CommandLineArguments::Keys::Required, L);
    cmd.readArray("n", CommandLineArguments::Keys::Required, n);
    d = L/n;
    bool isPBC = cmd.readArray("box", CommandLineArguments::Keys::Optional, box);

    // read the input file and set up the mesh 
    std::vector<Mesh> Meshs;
    for (auto fname : inputfname){
        Mesh m;
        MeshTools::readPLYlibr(fname, m);
        if (isPBC){
            m.setBoxLength(box); 
        }
        Meshs.push_back(m);
    }

    // initialize the points
    std::vector<Real3> points;

    // populate the points 
    for (int i=0;i<n[0];i++)
    {
        for (int j=0;j<n[1];j++)
        {
            Real3 p = {{height_pos, height_pos, height_pos}};
            p[ProjectedIndex[0]] = (i+0.5) * d[0];
            p[ProjectedIndex[1]] = (j+0.5) * d[1];

            points.push_back(p);
        }
    }

    // whether or not we hit the meshs
    std::vector<bool> hitRefMesh(points.size(), false);
    std::vector<Real3> hitNormal(points.size());

    // see if any points hits any of the reference mesh 
    for (auto& m : Meshs)
    {
        m.CalcVertexNormals();
        const auto& refFace = m.gettriangles();
        const auto& refverts= m.getvertices();
        const auto& refNormal=m.getFaceNormals();

        // fix the triangles
        std::vector<std::array<Real3,3>> fixed_tri_pbc(refFace.size());
        #pragma omp parallel for
        for (int i=0;i<refFace.size();i++)
        {
            std::array<Real3,3> t;
            Real3 A, B, C;
            A = refverts[refFace[i].triangleindices_[0]].position_;
            B = refverts[refFace[i].triangleindices_[1]].position_;
            C = refverts[refFace[i].triangleindices_[2]].position_;

            // if the mesh is periodic, then we shift the other 2 vertices with respect to the first vertex 
            if (isPBC)
            {
                MeshTools::ShiftPeriodicTriangle(const_cast<std::vector<vertex>&>(refverts), \
                                                const_cast<INT3&>(refFace[i].triangleindices_), box, A, B, C);
            }

            fixed_tri_pbc[i] = {{A,B,C}};
        }

        #pragma omp parallel for
        for (int i=0;i<points.size();i++)
        {
            Real3 O = points[i];
            for (int j=0;j<refFace.size();j++)
            {
                Real3 projectedD;
                std::array<Real3,3> tri = fixed_tri_pbc[j];
                Real3 A, B, C;
                Real t, u ,v;
                A = tri[0];
                B = tri[1];
                C = tri[2];

                bool isIntersect = MeshTools::MTRayTriangleIntersection(A, B, C, \
                                                    O, D, t, u ,v);

                if (isIntersect)
                {
                    hitRefMesh[i] = true;
                    hitNormal[i]  = refNormal[j];
                    break;
                }
            }
        }
    }

    // start to output file
    std::ofstream ofs(outputfname);
    int index=0;
    for (int i=0;i<n[0];i++)
    {
        for (int j=0;j<n[1];j++)
        {
            ofs << i << " " << j << " " << hitRefMesh[index] << "\n";
            index++;
        }
    }

    ofs.close();

    // if we are doing normal mapping 
    if (normalMap){
        std::ofstream ofs;
        ofs.open(normalMapfname);
        int index=0;
        for (int i=0;i<n[0];i++)
        {
            for (int j=0;j<n[1];j++)
            {
                Real3 normal = hitNormal[index];
                Real R = (normal[ProjectedIndex[0]] + 1)/2;
                Real G = (normal[ProjectedIndex[1]] + 1)/2;
                Real B = (normal[0] + 1)/2;

                ofs << i << " " << j << " " << R << " " << G << " " << B << "\n";

                index++;
            }
        }
        ofs.close();
    }
}

void MeshActions::Project3dCurvature(CommandLineArguments& cmd)
{
    using INT2 = std::array<int,2>;

    // define input output file names 
    std::string inputfname, outputfname, fcfname;
    std::vector<std::string> refMeshfname;

    // some large number 
    Real origin_pos=100.0;
    Real2 L, n, d, ProjectedIndex;
    std::vector<std::vector<Real>> fc;
    Real3 box, D;

    // read the inputs 
    cmd.readValue("i", CommandLineArguments::Keys::Required, inputfname);
    cmd.readValue("o", CommandLineArguments::Keys::Optional, outputfname);
    cmd.readArray("RayDirection", CommandLineArguments::Keys::Required,D);
    LinAlg3x3::normalize(D);
    cmd.readArray("ProjectedIndex", CommandLineArguments::Keys::Required, ProjectedIndex);

    // reference mesh
    cmd.readVector("refMesh", CommandLineArguments::Keys::Optional, refMeshfname);

    // the origin 
    cmd.readValue("origin", CommandLineArguments::Keys::Optional, origin_pos);
    cmd.readArray("L", CommandLineArguments::Keys::Required, L);
    cmd.readArray("n", CommandLineArguments::Keys::Required, n);
    d = L/n;
    cmd.readValue("FaceCurvatureFile", CommandLineArguments::Keys::Required, fcfname);
    bool isPBC = cmd.readArray("box", CommandLineArguments::Keys::Optional, box);

    // read the tabulated data for curvature
    StringTools::ReadTabulatedData(fcfname, fc);
    int colsize = fc[0].size();

    // read the input file and set up the mesh 
    Mesh mesh;
    MeshTools::readPLYlibr(inputfname, mesh);
    if (isPBC){mesh.setBoxLength(box);}

    // read the reference meshes
    std::vector<Mesh> refMesh;
    for (auto fname : refMeshfname){
        Mesh ref;
        MeshTools::readPLYlibr(fname, ref);
        if (isPBC){ref.setBoxLength(box); }
        refMesh.push_back(ref);
    }


    // initialize and populate the lattice points
    std::vector<Real3> lattice;
    for (int i=0;i<n[0];i++){
        for (int j=0;j<n[1];j++){
            Real3 p = {{origin_pos, origin_pos, origin_pos}};
            p[ProjectedIndex[0]] = (i+0.5) * d[0];
            p[ProjectedIndex[1]] = (j+0.5) * d[1];

            lattice.push_back(p);
        }
    }

    // points curvature
    std::vector<std::vector<Real>> pointCurvature(lattice.size());
    std::vector<bool> hitRefMesh(lattice.size(), false);
    std::vector<Real> hitRefMeshHeightMin(lattice.size(), 1e10);
    std::vector<std::vector<Real>> hitRefMeshHeightTotal(lattice.size());

    // see if any points hits any of the reference mesh 
    for (auto& m : refMesh){
        // for each reference mesh --> calculate vertex normal
        m.CalcVertexNormals();

        // get the properties of the mesh
        const auto& refFace = m.gettriangles();
        const auto& refverts= m.getvertices();
        const auto& refNormal=m.getFaceNormals();

        #pragma omp parallel for
        for (int i=0;i<lattice.size();i++){
            // define the origin points
            Real3 O = lattice[i];
            for (int j=0;j<refFace.size();j++){
                Real3 projectedD;
                Real3 A, B, C;
                A = refverts[refFace[j].triangleindices_[0]].position_;
                B = refverts[refFace[j].triangleindices_[1]].position_;
                C = refverts[refFace[j].triangleindices_[2]].position_;

                Real t, u, v;

                // if the mesh is periodic, then we shift the other 2 vertices with respect to the first vertex 
                if (m.isPeriodic()){
                    MeshTools::ShiftPeriodicTriangle(const_cast<std::vector<vertex>&>(refverts), \
                                                    const_cast<INT3&>(refFace[j].triangleindices_), box, A, B, C);
                }

                // hit ref mesh height
                if (MeshTools::MTRayTriangleIntersection(A, B, C, O, D, t, u, v)){
                    hitRefMeshHeightTotal[i].push_back(t);
                }
            }
        }
    }

    #pragma omp parallel for
    for (int i=0;i<lattice.size();i++){
        if (hitRefMeshHeightTotal[i].size() != 0){
            hitRefMeshHeightMin[i] = Algorithm::min(hitRefMeshHeightTotal[i]);
        }
    }

    // faces, verts of the mesh 
    mesh.CalcVertexNormals();
    const auto& face = mesh.gettriangles();
    const auto& verts= mesh.getvertices();
    const auto& faceNormal = mesh.getFaceNormals();

    // shift the face
    std::vector<std::array<Real3,3>> shiftedTriangle(face.size());

    #pragma omp parallel for 
    for (int i=0;i<face.size();i++)
    {
        Real3 A, B, C;
        // if the mesh is periodic, then we shift the other 2 vertices with respect to the first vertex 
        if (isPBC){
            MeshTools::ShiftPeriodicTriangle(const_cast<std::vector<vertex>&>(verts), \
                                            const_cast<INT3&>(face[i].triangleindices_), box, A, B, C);
        }
        else{
            A = verts[face[i].triangleindices_[0]].position_;
            B = verts[face[i].triangleindices_[1]].position_;
            C = verts[face[i].triangleindices_[2]].position_;
        }

        shiftedTriangle[i] = {{A,B,C}};
    }

    // iterate over the points 
    #pragma omp parallel for 
    for (int i=0;i<lattice.size();i++){
        Real3 O = lattice[i];
        std::vector<int> faceIndex;
        std::vector<Real> value;

        for (int j=0;j<face.size();j++){
            Real3 projectedD;
            Real t, u, v;
            Real3 A=shiftedTriangle[j][0], B = shiftedTriangle[j][1], C=shiftedTriangle[j][2];

            // shift centroid with respect to the O
            Real3 centroid = 1.0 / 3.0 * (A + B + C);
            Real3 shiftVec;
            mesh.CalculateShift(centroid, O, shiftVec);
            for (int i=0;i<2;i++){
                A[ProjectedIndex[i]] += shiftVec[ProjectedIndex[i]];
                B[ProjectedIndex[i]] += shiftVec[ProjectedIndex[i]];
                C[ProjectedIndex[i]] += shiftVec[ProjectedIndex[i]];
            }

            if (MeshTools::MTRayTriangleIntersection(A,B,C,O,D,t,u,v)){
                faceIndex.push_back(j);
                value.push_back(t);
            }
        }

        // value size = 0 means that this particular Ray did not hit any of the triangles
        std::vector<Real> temp(colsize,-100);
        if (value.size() != 0){
            // find the min element index of t --> the shorted the t, the closer it is to the viewing point
            int minIndex = Algorithm::argmin(value);
            if (hitRefMeshHeightMin[i] != 1e10){
                pointCurvature[i] = temp;
            }
            else{
                int fcIndex = faceIndex[minIndex];
                pointCurvature[i] = fc[fcIndex];
            }
            // if (value[minIndex] < hitRefMeshHeightMin[i]){
            //     int fcIndex = faceIndex[minIndex];
            //     pointCurvature[i] = fc[fcIndex];
            // }
            // else {
            //     pointCurvature[i] = temp;
            // }
        }
        else{
            pointCurvature[i] = temp;
        }
    }

    // start to output file
    std::ofstream ofs(outputfname);
    int index=0;
    for (int i=0;i<n[0];i++){
        for (int j=0;j<n[1];j++){
            ofs << i << " " << j << " " << pointCurvature[index] << "\n";
            index++;
        }
    }

    ofs.close();
}

void MeshActions::MeshDistanceCutoff(CommandLineArguments& cmd)
{
    std::string inputfname, reffname, indfname="index.out";
    std::string outputfname="cut.ply";

    // distance of hydration shell in [nm]
    Real cutoff=0.3, cutoffsq;
    Real3 box;

    // distance index
    std::vector<int> distanceIndex = {{0,1,2}};

    bool writePBC=true;

    // read the inputs 
    cmd.readValue("i", CommandLineArguments::Keys::Required, inputfname);
    cmd.readValue("o", CommandLineArguments::Keys::Optional, outputfname);
    cmd.readValue("ref", CommandLineArguments::Keys::Required, reffname);
    cmd.readValue("cutoff", CommandLineArguments::Keys::Required, cutoff);
    cmd.readVector("distanceIndex", CommandLineArguments::Keys::Optional, distanceIndex);
    cmd.readBool("writePBC", CommandLineArguments::Keys::Optional, writePBC);
    bool indread = cmd.readValue("ind", CommandLineArguments::Keys::Optional, indfname);
    cutoffsq = cutoff * cutoff;
    bool isPBC = cmd.readArray("box", CommandLineArguments::Keys::Optional,box);

    // set up the mesh object 
    Mesh mesh, refmesh;
    MeshTools::readPLYlibr(inputfname, mesh);
    MeshTools::readPLYlibr(reffname, refmesh);

    // calculate vertices differences 
    const auto& v = mesh.getvertices();
    const auto& f = mesh.gettriangles();
    const auto& refv = refmesh.getvertices();
    std::vector<bool> indicator(v.size(), true);

    #pragma omp parallel for
    for (int i=0;i<v.size();i++)
    {
        for (int j=0;j<refv.size();j++)
        {
            Real3 distance;

            if (isPBC)
            {
                Real d;
                MeshTools::calculateDistance(v[i].position_, refv[j].position_, box, distance, d);
            }
            else
            {
                for (int k=0;k<3;k++)
                {
                    distance[k] = v[i].position_[k] - refv[j].position_[k];
                }
            }

            // calculate the squared distances 
            Real distsq = 0.0;
            for (int k=0;k<distanceIndex.size();k++)
            {
                distsq += distance[distanceIndex[k]] * distance[distanceIndex[k]];
            }

            if (distsq < cutoffsq)
            {
                indicator[i] = false;
                break;
            }
        }
    }

    std::vector<Real3> newV;
    std::vector<INT3> newF;
    std::vector<int> NewTIndices;
    std::vector<int> MapOldIndicesToNew(v.size(),-1);
    int index=0;

    // make new vertices and convert the old indices into new 
    for (int i=0;i<v.size();i++)
    {
        if (indicator[i])
        {
            newV.push_back(v[i].position_);
            MapOldIndicesToNew[i] = index;
            index++;
        }
    }

    // make new faces 
    for (int i=0;i<f.size();i++)
    {
        INT3 ind = f[i].triangleindices_;
        bool validTriangle=true;
        for (int j=0;j<3;j++)
        {
            if (! indicator[ind[j]])
            {
                validTriangle=false;
            }
        }

        if (validTriangle)
        {
            INT3 newT;
            for (int j=0;j<3;j++)
            {
                newT[j] = MapOldIndicesToNew[ind[j]];
            }
            
            if (! writePBC)
            {
                if (! MeshTools::IsPeriodicTriangle(const_cast<std::vector<vertex>&>(v), const_cast<INT3&>(f[i].triangleindices_), box))
                {
                    NewTIndices.push_back(i);
                    newF.push_back(newT);
                }
            }
            else
            {
                NewTIndices.push_back(i);
                newF.push_back(newT);
            }
        }
    }

    // output the ply file
    MeshTools::writePLY(outputfname, newV, newF);

    // output the indices 
    if (indread)
    {
        std::ofstream ofs;
        ofs.open(indfname);
        for (int i=0;i<NewTIndices.size();i++)
        {
            ofs << NewTIndices[i] << " ";
        }
        ofs.close();
    }
}

void MeshActions::MinimumMeshDistance(CommandLineArguments& cmd)
{
    std::string inputfname, outputfname, reffname;
    Real3 box;
    std::vector<int> DistanceIndices={{0,1,2}};

    // read in the inputs 
    cmd.readValue("i", CommandLineArguments::Keys::Required, inputfname);
    cmd.readValue("o", CommandLineArguments::Keys::Optional, outputfname);
    cmd.readValue("ref", CommandLineArguments::Keys::Required, reffname);
    cmd.readVector("index", CommandLineArguments::Keys::Optional, DistanceIndices);
    bool isPBC = cmd.readArray("box", CommandLineArguments::Keys::Optional, box);

    // mesh object 
    Mesh mesh, refmesh;
    MeshTools::readPLYlibr(inputfname, mesh);
    MeshTools::readPLYlibr(reffname, refmesh);
    if (isPBC){
        mesh.setBoxLength(box);refmesh.setBoxLength(box);
    }

    const auto& v = mesh.getvertices();
    const auto& refv = refmesh.getvertices();

    std::vector<Real> minDist(v.size(), 0.0);

    #pragma omp parallel for
    for (int i=0;i<v.size();i++)
    {
        std::vector<Real> distance(refv.size(),0.0);
        for (int j=0;j<refv.size();j++){
            Real distsq=0.0;
            if (isPBC){
                Real3 r;
                MeshTools::calculateDistance(v[i].position_, refv[j].position_, box, r, distsq);
                distsq = 0.0;

                for (int k=0;k<DistanceIndices.size();k++){
                    distsq += r[DistanceIndices[k]] * r[DistanceIndices[k]];
                }
            }
            else{
                for (int k=0;k<DistanceIndices.size();k++){
                    distsq += std::pow((v[i].position_[DistanceIndices[k]]-refv[j].position_[DistanceIndices[k]]),2.0);
                }
            }
            Real dist = std::sqrt(distsq);
            distance[j] = dist;
        }
        Real min_dist = *std::min_element(distance.begin(), distance.end());
        minDist[i] = min_dist;
    }

    std::ofstream ofs;
    ofs.open(outputfname);
    for (int i=0;i<v.size();i++){
        ofs << i << " " << minDist[i] << "\n";
    }
    ofs.close();
}

void MeshActions::ClipMesh(CommandLineArguments& cmd){
    using Intersector = MeshPlaneIntersect<Real, int>;

    std::vector<Intersector::Vec3D> vertices, normals;
    std::vector<Intersector::Face> faces;
    Real3 p, o;
    int skip=1;

    std::string inputfname;
    std::string outputfname="intersect.out";

    cmd.readValue("o", CommandLineArguments::Keys::Optional, outputfname);
    cmd.readArray("plane", CommandLineArguments::Keys::Required, p);
    cmd.readArray("points", CommandLineArguments::Keys::Required, o);
    cmd.readValue("skip", CommandLineArguments::Keys::Optional, skip);

    Mesh mesh;
    MeshTools::readPLYlibr(inputfname, mesh);
    mesh.CalcVertexNormals();

    const auto& v = mesh.getvertices();
    const auto& f = mesh.gettriangles();

    // copy the data to mesh plane intersection code 
    for (int i=0;i<v.size();i++){
        vertices.push_back(v[i].position_);
        normals.push_back(v[i].normals_);
    }

    for (int i=0;i<f.size();i++){
        faces.push_back(f[i].triangleindices_);
    }

    // create the mesh for mesh-plane intersection
    Intersector::Mesh m(vertices, faces, normals);
    Intersector::Plane plane;
    plane.normal = p;
    plane.origin = o;

    auto result = m.Clip(plane);
    std::cout << "result size = " << result.size() << "\n";
    std::cout << "result[0] points size = " << result[0].points.size() << "\n";
}

void MeshActions::MeshPlaneIntersection(CommandLineArguments& cmd)
{
    using Intersector = MeshPlaneIntersect<Real, int>;

    std::vector<Intersector::Vec3D> vertices, normals;
    std::vector<Intersector::Face> faces;
    Real3 p, o;
    int skip=1;

    std::string inputfname;
    std::string outputfname="intersect.out";

    cmd.readValue("i", CommandLineArguments::Keys::Required, inputfname);
    cmd.readValue("o", CommandLineArguments::Keys::Optional, outputfname);
    cmd.readArray("plane", CommandLineArguments::Keys::Required, p);
    cmd.readArray("points", CommandLineArguments::Keys::Required, o);
    cmd.readValue("skip", CommandLineArguments::Keys::Optional, skip);

    Mesh mesh;
    MeshTools::readPLYlibr(inputfname, mesh);
    mesh.CalcVertexNormals();

    const auto& v = mesh.getvertices();
    const auto& f = mesh.gettriangles();

    // copy the data to mesh plane intersection code 
    for (int i=0;i<v.size();i++){
        vertices.push_back(v[i].position_);
        normals.push_back(v[i].normals_);
    }

    for (int i=0;i<f.size();i++){
        faces.push_back(f[i].triangleindices_);
    }

    // create the mesh for mesh-plane intersection
    Intersector::Mesh m(vertices, faces, normals);
    Intersector::Plane plane;
    plane.normal = p;
    plane.origin = o;

    auto result = m.Intersect(plane);
    if (result.size() < 1){
        std::cout << "No path was found." << "\n";
        return;
    }

    std::ofstream ofs;
    ofs.open(outputfname);
    ofs << "# vx vy vz nx ny nz\n";
    for (int i=0;i<result[0].points.size();i+=skip){
        for (int j=0;j<3;j++){
            ofs << result[0].points[i][j] << " ";
        }

        for (int j=0;j<3;j++){
            ofs << result[0].normals[i][j] << " ";
        }
        ofs << "\n";
    }
    ofs.close();
}

void MeshActions::FindVertexNeighbors(CommandLineArguments& cmd)
{
    std::string inputfname, outputfname="neighbors.out";

    int n = 1;
    cmd.readValue("i", CommandLineArguments::Keys::Required, inputfname);
    cmd.readValue("n", CommandLineArguments::Keys::Optional, n);
    cmd.readValue("o", CommandLineArguments::Keys::Optional, outputfname);

    // read in the mesh
    Mesh mesh;
    MeshTools::readPLYlibr(inputfname, mesh);

    // calculate the neighbor indices 
    std::vector<std::vector<int>> neighborInd;
    std::vector<std::vector<int>> neighborInd_nVertex;
    MeshTools::CalculateVertexNeighbors(mesh, neighborInd);
    Graph::getNearbyIndicesNVertexAway(neighborInd, n, neighborInd_nVertex);

    std::ofstream ofs;
    ofs.open(outputfname);

    for (int i=0;i<neighborInd.size();i++){
        ofs << i << " ";
        for (auto ind : neighborInd_nVertex[i]){
            ofs << ind << " ";
        }
        ofs << "\n";
    }

    ofs.close();
}

void MeshActions::FindBoundaryVertices(CommandLineArguments& cmd){
    std::string inputfname, outputfname="boundary.out";

    cmd.readValue("i", CommandLineArguments::Keys::Required, inputfname);
    cmd.readValue("o", CommandLineArguments::Keys::Optional, outputfname);

    Mesh mesh;
    MeshTools::readPLYlibr(inputfname, mesh);

    std::map<INT2, std::vector<int>> MapEdgeToFace;
    std::vector<std::vector<INT2>> MapVertexToEdge;
    std::vector<bool> BoundaryIndicator;

    MeshTools::MapEdgeToFace(mesh,MapEdgeToFace, MapVertexToEdge);
    MeshTools::CalculateBoundaryVertices(mesh, MapEdgeToFace,BoundaryIndicator);

    std::ofstream ofs;
    ofs.open(outputfname);

    for (int i=0;i<BoundaryIndicator.size();i++){
        if (BoundaryIndicator[i]){
            ofs << i << " ";
        }
    }

    ofs.close();
}

void MeshActions::CutOverlappedRegion(CommandLineArguments& cmd)
{
    using INT2 = std::array<int,2>;

    // define input output file names 
    std::string inputfname, reffname, normalMapfname, outputfname;

    // some large number 
    Real origin_pos=100.0;
    int colnum;
    Real2 L, n, d, ProjectedIndex;
    std::vector<Real> fc;
    Real3 box, D;

    // read the inputs 
    cmd.readValue("i", CommandLineArguments::Keys::Required, inputfname);
    cmd.readValue("ref", CommandLineArguments::Keys::Required, reffname);
    cmd.readValue("o", CommandLineArguments::Keys::Optional, outputfname);
    cmd.readArray("RayDirection", CommandLineArguments::Keys::Required,D);
    LinAlg3x3::normalize(D);
    bool isPBC = cmd.readArray("box", CommandLineArguments::Keys::Optional, box);

    // read the input file and set up the mesh 
    Mesh m, refm;
    MeshTools::readPLYlibr(inputfname, m);
    MeshTools::readPLYlibr(reffname, refm);
    if (isPBC){
        m.setBoxLength(box);
        refm.setBoxLength(box);
    }

    // see if any points hits any of the reference mesh 
    const auto& Face = m.gettriangles();
    const auto& verts= m.getvertices();
    const auto& refFace = refm.gettriangles();
    const auto& refverts= refm.getvertices();

    // fix the triangles
    std::vector<std::array<Real3,3>> fixed_tri_pbc(Face.size());
    std::vector<std::array<Real3,3>> fixed_reftri_pbc(refFace.size());

    #pragma omp parallel for
    for (int i=0;i<Face.size();i++)
    {
        std::array<Real3,3> t;
        Real3 A, B, C;
        A = verts[Face[i].triangleindices_[0]].position_;
        B = verts[Face[i].triangleindices_[1]].position_;
        C = verts[Face[i].triangleindices_[2]].position_;

        // if the mesh is periodic, then we shift the other 2 vertices with respect to the first vertex 
        if (isPBC){
            MeshTools::ShiftPeriodicTriangle(const_cast<std::vector<vertex>&>(verts), \
                                            const_cast<INT3&>(Face[i].triangleindices_), box, A, B, C);
        }

        fixed_tri_pbc[i] = {{A,B,C}};
    }

    #pragma omp parallel for
    for (int i=0;i<refFace.size();i++)
    {
        std::array<Real3,3> t;
        Real3 A, B, C;
        A = refverts[refFace[i].triangleindices_[0]].position_;
        B = refverts[refFace[i].triangleindices_[1]].position_;
        C = refverts[refFace[i].triangleindices_[2]].position_;

        // if the mesh is periodic, then we shift the other 2 vertices with respect to the first vertex 
        if (isPBC){
            MeshTools::ShiftPeriodicTriangle(const_cast<std::vector<vertex>&>(refverts), \
                                            const_cast<INT3&>(refFace[i].triangleindices_), box, A, B, C);
        }

        fixed_reftri_pbc[i] = {{A,B,C}};
    }

    std::vector<bool> keeptri(fixed_tri_pbc.size(), true);
    #pragma omp parallel for
    for (int i=0;i<fixed_tri_pbc.size();i++)
    {
        Real3 O = (fixed_tri_pbc[i][0] + fixed_tri_pbc[i][1] + fixed_tri_pbc[i][2]) * 1.0/3.0;
        for (int j=0;j<fixed_reftri_pbc.size();j++)
        {
            Real3 A,B,C;
            A = fixed_reftri_pbc[j][0], B=fixed_reftri_pbc[j][1], C=fixed_reftri_pbc[j][2];
            Real t, u ,v;

            if (MeshTools::MTRayTriangleIntersection(A,B,C,O,D,t,u,v)){
                keeptri[i] = false;
                break;
            }
        }
    }

    Mesh newMesh;
    auto& newF = newMesh.accesstriangles();
    auto& newV = newMesh.accessvertices();

    for (int i=0;i<Face.size();i++){
        if (keeptri[i]){
            newF.push_back(Face[i]);
        }
    }

    for (int i=0;i<verts.size();i++){
        newV.push_back(verts[i]);
    }
    std::vector<std::vector<int>> neighborIndex;
    MeshTools::CalculateVertexNeighbors(newMesh, neighborIndex);

    int index=0;
    std::vector<int> MapOldToNewIndex(newV.size(),-1);
    for (int i=0;i<neighborIndex.size();i++){
        if (neighborIndex[i].size() != 0){
            MapOldToNewIndex[i] = index;            
            index++;
        }
    }

    std::vector<INT3> newFace;
    std::vector<Real3> newVerts;
    for (int i=0;i<newF.size();i++){
        INT3 ind;
        for (int j=0;j<3;j++){
            ind[j] = MapOldToNewIndex[newF[i].triangleindices_[j]];
        }
        newFace.push_back(ind);
    }

    for (int i=0;i<newV.size();i++){
        if (neighborIndex[i].size() != 0){
            newVerts.push_back(verts[i].position_);
        }
    }

    MeshTools::writePLY(outputfname, newVerts, newFace);
}

void MeshActions::DecimateDegenerateTriangles(CommandLineArguments& cmd)
{
    std::string inputfname, outputfname="out.ply";
    Real3 box;

    cmd.readValue("i", CommandLineArguments::Keys::Required, inputfname);
    cmd.readValue("o", CommandLineArguments::Keys::Optional, outputfname);
    bool isPBC = cmd.readArray("box", CommandLineArguments::Keys::Optional, box);

    Mesh m;
    MeshTools::readPLYlibr(inputfname, m);

    if (isPBC){m.setBoxLength(box);}

    MeshTools::decimateDegenerateTriangle(m);

    MeshTools::writePLY(outputfname, m);
}

void MeshActions::RefineBoundary(CommandLineArguments& cmd){
    std::string inputfname, outputfname="evolved.ply", FixedIndexFile, neighbors, CurvatureStepsize, k0, maxCurvaturestep="1e5", tolerance, surface_shape_name, boundaryZ_offset;
    std::unique_ptr<AFP_shape> shape;
    Real3 box;
    Real contact_angle_val=-0.5;
    Real BoundaryStepSize=0.01;

    // initialize evolution
    refineptr r;
    ParameterPack EvolutionPack, curvePack, shapePack;
    bool fairing=false, cleanMesh=false;
    std::string curve_type="curvefit",debug="false";

    cmd.readString("i", CommandLineArguments::Keys::Required, inputfname);
    cmd.readString("o", CommandLineArguments::Keys::Optional, outputfname);

    // use curve fit 
    bool isPBC = cmd.readArray("box", CommandLineArguments::Keys::Optional, box);
    if (curve_type == "curvefit"){
        cmd.readString("neighbors", CommandLineArguments::Keys::Required, neighbors);
        // fill parameter pack
        curvePack.insert("neighbors", neighbors);
        curvePack.insert("type", "curvefit");
        curvePack.insert("name", "temp");
    }
    else if (curve_type == "jetfit"){
        std::string MongeCoefficient, degree;
        cmd.readString("neighbors", CommandLineArguments::Keys::Required, neighbors);
        cmd.readString("MongeCoefficient", CommandLineArguments::Keys::Required, MongeCoefficient);
        cmd.readString("degree", CommandLineArguments::Keys::Required, degree);

        curvePack.insert("neighbors", neighbors);
        curvePack.insert("type", "jetfit");
        curvePack.insert("name", "temp");
        curvePack.insert("MongeCoefficient", MongeCoefficient);
        curvePack.insert("degree", degree);
    }
    else{
        ASSERT((true==false), "The curvature type " << curve_type << " is not valid.");
    }

    // start reading the step size information etc.
    cmd.readString("CurvatureStep", CommandLineArguments::Keys::Required, CurvatureStepsize);
    cmd.readString("k0", CommandLineArguments::Keys::Required, k0);
    cmd.readString("maxCurvaturestep", CommandLineArguments::Keys::Optional, maxCurvaturestep);
    cmd.readString("tolerance", CommandLineArguments::Keys::Optional, tolerance);
    cmd.readBool("fairing",CommandLineArguments::Keys::Optional, fairing);
    cmd.readBool("cleanMesh", CommandLineArguments::Keys::Optional, cleanMesh);

    // insert everything into evolution pack.
    EvolutionPack.insert("Curvature", curvePack);
    EvolutionPack.insert("name", "refine");
    EvolutionPack.insert("stepsize", CurvatureStepsize);
    EvolutionPack.insert("k0", k0);
    EvolutionPack.insert("maxstep", maxCurvaturestep);
    EvolutionPack.insert("tolerance", tolerance);
    EvolutionPack.insert("debug", debug);

    // read the fixed index file
    if (cmd.readString("fixed_index_file", CommandLineArguments::Keys::Optional, FixedIndexFile)){
        EvolutionPack.insert("fixed_index_file", FixedIndexFile);
    }
        
    if (cleanMesh){
        std::string edgelength;
        cmd.readString("edgeLengthCutoff", CommandLineArguments::Keys::Required, edgelength);
        EvolutionPack.insert("cleanMesh", "true");
        EvolutionPack.insert("edgeLengthCutoff", edgelength);
    }

    if (fairing){
        std::string fairing_iteration, fairing_step;
        cmd.readString("fairing_iteration", CommandLineArguments::Keys::Required, fairing_iteration);
        cmd.readString("fairing_step", CommandLineArguments::Keys::Required, fairing_step);
        EvolutionPack.insert("fairing", "true");
        EvolutionPack.insert("fairing_iteration", fairing_iteration);
        EvolutionPack.insert("fairing_step", fairing_iteration);
    }

    // initialize shape pack

    cmd.readString("shape", CommandLineArguments::Keys::Required, surface_shape_name);
    if (surface_shape_name == "superegg"){
        std::string a,b,n,zmax,a_taper,b_taper,a_alpha,b_alpha;
        std::vector<std::string> center;
        bool read;
        cmd.readValue("a", CommandLineArguments::Keys::Required, a);
        shapePack.insert("a", a);
        cmd.readValue("b", CommandLineArguments::Keys::Required, b);
        shapePack.insert("b", b);
        cmd.readValue("zmax", CommandLineArguments::Keys::Required, zmax);
        shapePack.insert("zmax", zmax);
        cmd.readVector("center", CommandLineArguments::Keys::Required, center);
        shapePack.insert("center", center);
        read = cmd.readValue("n", CommandLineArguments::Keys::Optional, n);
        if (read){
            shapePack.insert("n", n);
        }
        read = cmd.readValue("boundaryZ_offset", CommandLineArguments::Keys::Optional, boundaryZ_offset);
        if (read){
            shapePack.insert("offset_height", boundaryZ_offset);
        }
        read = cmd.readValue("a_taper", CommandLineArguments::Keys::Optional, a_taper);
        if (read){
            shapePack.insert("a_taper", a_taper);
        }
        read = cmd.readValue("b_taper", CommandLineArguments::Keys::Optional, b_taper);
        if (read){
            shapePack.insert("b_taper", b_taper);
        }
        read = cmd.readValue("a_alpha", CommandLineArguments::Keys::Optional, a_alpha);
        if (read){
            shapePack.insert("a_alpha", a_alpha);
        }
        read = cmd.readValue("b_alpha", CommandLineArguments::Keys::Optional, b_alpha);
        if (read){
            shapePack.insert("b_alpha", b_alpha);
        }

        // initialize the shape 
        shape = std::make_unique<SuperEgg>(shapePack);
    }
    else{
        ASSERT((true == false), "The shape " << surface_shape_name << " is not implemented yet.");
    }


    // read boundary refinement information
    int maxBoundaryStep, maxKStep=100;
    Real contact_angle = -0.5, std_limit=1e-3;
    Real k_step=0.1, k;
    Real k_tolerance=0.001;
    k = StringTools::StringToType<Real>(k0);
    std::string mean_angle_output, kappa_output;
    bool output_mean_angle=false, output_kappa=false;

    cmd.readValue("BoundaryStep", CommandLineArguments::Keys::Optional, BoundaryStepSize);
    cmd.readValue("maxBoundaryStep", CommandLineArguments::Keys::Optional, maxBoundaryStep);
    cmd.readValue("ContactAngle", CommandLineArguments::Keys::Optional, contact_angle);
    cmd.readValue("std_limit", CommandLineArguments::Keys::Optional, std_limit);
    cmd.readValue("k_step", CommandLineArguments::Keys::Optional, k_step);
    cmd.readValue("maxk_step", CommandLineArguments::Keys::Optional, maxKStep);
    cmd.readValue("k_tolerance", CommandLineArguments::Keys::Optional, k_tolerance);
    output_mean_angle = cmd.readValue("mean_angle_output", CommandLineArguments::Keys::Optional, mean_angle_output);
    output_kappa = cmd.readValue("kappa_output", CommandLineArguments::Keys::Optional, kappa_output);

    // refine pointer 
    MeshRefineStrategyInput input = {{EvolutionPack}};
    r = refineptr(MeshRefineStrategyFactory::Factory::instance().create("CurvatureEvolution", input));


    // read input mesh 
    Mesh m, original_m;
    MeshTools::readPLYlibr(inputfname, m);
    original_m = m;
    if (isPBC){m.setBoxLength(box);original_m.setBoxLength(box);}

    // boundary information
    std::vector<bool> boundaryIndicator;
    MeshTools::CalculateBoundaryVertices(m, boundaryIndicator);

    // first get the curvature to about the value predicted to be -0.5
    // auto rr = dynamic_cast<CurvatureEvolution*>(r.get());
    // ASSERT((rr!= nullptr), "Conversion wrong");
    // rr->setMeanCurvature(k);
    // rr->refine(original_m);

    std::vector<Real> cangle_list;
    bool exceeded_shape=false;

    // start iterating
    for (int i=0;i<maxKStep;i++){
        // set m back to the original mesh
        m = original_m;

        // set the curvature to the curvature value
        auto rr = dynamic_cast<CurvatureEvolution*>(r.get());
        ASSERT((rr!= nullptr), "Conversion wrong");
        rr->setMeanCurvature(k);
        std::cout << "We are at curvature " << k << std::endl;

        // obtain the vertices 
        auto& vertices = m.accessvertices();

        // whether or not it is converged
        bool converged = false;

        // refine it first such that it is at that curvature
        rr->refine(m);

        Real curr_mean = 1e10;
        Real curr_std  = 1e10;
        Real mean;

        for (int j=0;j<maxBoundaryStep;j++){
            cangle_list.clear();
            std::vector<Real3> s_vec;
            std::vector<Real3> n_vec;
            for (int k=0;k<boundaryIndicator.size();k++){
                if (MeshTools::IsBoundary(k, boundaryIndicator)){
                    double3 v, up_v, down_v;
                    Real3 s, t, up_t, up_s, down_t, down_s;
                    v[0] = vertices[k].position_[0];
                    v[1] = vertices[k].position_[1];
                    v[2] = vertices[k].position_[2];

                    // check if we can calculate the normal and tangent
                    bool b =shape->CalculateNormalAndTangent(v, t, s); 

                    // if not, then we set the exceed_shape to be true
                    if (! b){
                        exceeded_shape = true;
                        break;
                    }

                    Real cangle = LinAlg3x3::DotProduct(s, vertices[k].normals_);
                    Real diff   = cangle - contact_angle;
                    cangle_list.push_back(cangle);
                    s_vec.push_back(s);
                    n_vec.push_back(vertices[k].normals_);

                    // update the vertices by the tangent 
                    // move up in t
                    vertices[k].position_ = vertices[k].position_ - BoundaryStepSize * diff * t;
                }
            }

            if (exceeded_shape){
                std::cout << "Failed to find the solution." << std::endl;
                break;
            }

            Real var = Algorithm::calculateVariance(cangle_list);
            mean= Algorithm::calculateMean(cangle_list);
            Real diff = std::abs(mean - contact_angle);


            std::cout << "std = " << std::sqrt(var) << std::endl;
            std::cout << "mean = " << mean << std::endl;

            // break if the var is small and mean is small 
            if (std::sqrt(var) < std_limit){
                if (std::abs(mean - contact_angle) < k_tolerance){
                    std::cout << "std " << std::sqrt(var) << " is lower than the set value of " << std_limit << " and mean has converged and is " << mean << std::endl;
                    converged = true;
                }
                else{
                    std::cout << "std " << std::sqrt(var) << " is lower than the set value of " << std_limit << " but mean has not converged and is " << mean << std::endl;
                }

                break;
            }

            // update the current mean
            // if (std::sqrt(var) > curr_std){
            //    std::cout << "std starts to increase" << std::endl;
            //    break;
            // }

            curr_std = std::sqrt(var);

            // refine again after updating the boundary
            rr->refine(m);
        }

        if (exceeded_shape){
            break;
        }

        // quit iteration if converged = true
        if (converged == true){
            std::cout << "Converged." << std::endl;
            break;
        }

        std::cout << "mean - contact angle = " << mean - contact_angle << std::endl;
        k +=  (mean - contact_angle) * k_step;
    }

    // if curvature, write curvature
    if (output_kappa){
        std::ofstream ofs;
        ofs.open(kappa_output);
        ofs << k << "\n";
        ofs.close();
    }

    // if mean_angle, output mean angle
    const auto& v = m.getvertices();
    if (output_mean_angle){
        std::ofstream ofs;
        ofs.open(mean_angle_output);
        int ind = 0;
        for (int i=0;i<boundaryIndicator.size();i++){
            if (MeshTools::IsBoundary(i, boundaryIndicator)){
                ofs << v[i].position_[0] << " " << v[i].position_[1] << " " << v[i].position_[2] << " " << cangle_list[ind] << "\n"; 
                ind++;
            }
        }
        ofs.close();
    }
    

    MeshTools::writePLY(outputfname, m);
}

void MeshActions::CurvatureEvolution1(CommandLineArguments& cmd)
{
    std::string inputfname, outputfname="evolved.ply", FixedIndexFile, neighbors, stepsize, k0, maxstep="1e5", tolerance;
    Real3 box;

    // initialize curve ptr
    refineptr r;
    ParameterPack EvolutionPack, curvePack;
    bool fairing=false, cleanMesh=false;
    std::string curve_type="curvefit",debug="false";

    // read input output
    cmd.readString("i", CommandLineArguments::Keys::Required, inputfname);
    cmd.readString("o", CommandLineArguments::Keys::Optional, outputfname);
    cmd.readString("curve_type", CommandLineArguments::Keys::Optional, curve_type);
    cmd.readString("debug", CommandLineArguments::Keys::Optional, debug);

    // use curve fit 
    bool isPBC = cmd.readArray("box", CommandLineArguments::Keys::Optional, box);
    if (curve_type == "curvefit"){
        cmd.readString("neighbors", CommandLineArguments::Keys::Required, neighbors);
        // fill parameter pack
        curvePack.insert("neighbors", neighbors);
        curvePack.insert("type", "curvefit");
        curvePack.insert("name", "temp");
    }
    else if (curve_type == "jetfit"){
        std::string MongeCoefficient, degree;
        cmd.readString("neighbors", CommandLineArguments::Keys::Required, neighbors);
        cmd.readString("MongeCoefficient", CommandLineArguments::Keys::Required, MongeCoefficient);
        cmd.readString("degree", CommandLineArguments::Keys::Required, degree);

        curvePack.insert("neighbors", neighbors);
        curvePack.insert("type", "jetfit");
        curvePack.insert("name", "temp");
        curvePack.insert("MongeCoefficient", MongeCoefficient);
        curvePack.insert("degree", degree);
    }
    else{
        ASSERT((true==false), "The curvature type " << curve_type << " is not valid.");
    }

    // start reading the step size information etc.
    cmd.readString("stepsize", CommandLineArguments::Keys::Required, stepsize);
    cmd.readString("k0", CommandLineArguments::Keys::Required, k0);
    cmd.readString("maxstep", CommandLineArguments::Keys::Optional, maxstep);
    cmd.readString("tolerance", CommandLineArguments::Keys::Optional, tolerance);
    cmd.readBool("fairing",CommandLineArguments::Keys::Optional, fairing);
    cmd.readBool("cleanMesh", CommandLineArguments::Keys::Optional, cleanMesh);

    // insert everything into evolution pack.
    EvolutionPack.insert("Curvature", curvePack);
    EvolutionPack.insert("name", "refine");
    EvolutionPack.insert("stepsize", stepsize);
    EvolutionPack.insert("k0", k0);
    EvolutionPack.insert("maxstep", maxstep);
    EvolutionPack.insert("tolerance", tolerance);
    EvolutionPack.insert("debug", debug);

    // read the fixed index file
    if (cmd.readString("fixed_index_file", CommandLineArguments::Keys::Optional, FixedIndexFile)){
        EvolutionPack.insert("fixed_index_file", FixedIndexFile);
    }
        
    if (cleanMesh){
        std::string edgelength;
        cmd.readString("edgeLengthCutoff", CommandLineArguments::Keys::Required, edgelength);
        EvolutionPack.insert("cleanMesh", "true");
        EvolutionPack.insert("edgeLengthCutoff", edgelength);
    }

    if (fairing){
        std::string fairing_iteration, fairing_step;
        cmd.readString("fairing_iteration", CommandLineArguments::Keys::Required, fairing_iteration);
        cmd.readString("fairing_step", CommandLineArguments::Keys::Required, fairing_step);
        EvolutionPack.insert("fairing", "true");
        EvolutionPack.insert("fairing_iteration", fairing_iteration);
        EvolutionPack.insert("fairing_step", fairing_iteration);
    }

    // refine pointer 
    MeshRefineStrategyInput input = {{EvolutionPack}};
    r = refineptr(MeshRefineStrategyFactory::Factory::instance().create("CurvatureEvolution", input));

    // read input mesh 
    Mesh m;
    // set box length for mesh
    MeshTools::readPLYlibr(inputfname, m);
    if (isPBC){m.setBoxLength(box);}

    // refine the m1.0
    r->refine(m);

    // output the mesh
    MeshTools::writePLY(outputfname, m);
}

void MeshActions::FindIsolatedFace(CommandLineArguments& cmd){
    std::string inputfname, outputfname="isolated.out";

    cmd.readString("i", CommandLineArguments::Keys::Required, inputfname);
    cmd.readString("o", CommandLineArguments::Keys::Optional, outputfname);

    Mesh m;
    MeshTools::readPLYlibr(inputfname, m);

    const auto& t = m.gettriangles();
    int nf = t.size();
    std::vector<INT3> IsolatedTriangles;
    std::map<INT2, std::vector<int>> mapEdgeToFace;
    std::vector<std::vector<INT2>> mapVertexToEdge;
    MeshTools::MapEdgeToFace(m, mapEdgeToFace, mapVertexToEdge);

    for (int i=0;i<nf;i++){
        if (MeshTools::IsIsolatedFace(m, i, mapEdgeToFace)){
            IsolatedTriangles.push_back(t[i].triangleindices_);
        }
    }

    std::ofstream ofs;
    ofs.open(outputfname);
    for (int i=0;i<IsolatedTriangles.size();i++){
        ofs << IsolatedTriangles[i] << "\n";
    }

    ofs.close();
}

void MeshActions::CutTeethlikeFace(CommandLineArguments& cmd){
    std::string inputfname, outputfname="cut.ply";

    cmd.readString("i", CommandLineArguments::Keys::Required, inputfname);
    cmd.readString("o", CommandLineArguments::Keys::Optional, outputfname);

    Mesh m;
    MeshTools::readPLYlibr(inputfname, m);

    bool teethlike = true;
    int max = 5;
    int iter=1;

    while (teethlike){
        const auto& t = m.gettriangles();
        int nf = t.size();

        std::vector<int> teethlikeTriangles;
        std::map<INT2, std::vector<int>> mapEdgeToFace;
        std::vector<std::vector<INT2>> mapVertexToEdge;
        std::vector<triangle> newT;
        MeshTools::MapEdgeToFace(m, mapEdgeToFace, mapVertexToEdge);

        for (int i=0;i<nf;i++){
            if (MeshTools::IsTeethlikeFace(m, i, mapEdgeToFace)){
                teethlikeTriangles.push_back(i);
            }
            else{
                newT.push_back(t[i]);
            }
        }

        if (teethlikeTriangles.size() == 0){
            teethlike = false;
        }
        else{
            auto& oldT = m.accesstriangles();
            oldT.clear();
            oldT.insert(oldT.end(), newT.begin(), newT.end());
        }

        iter++;
        if (iter > max){
            std::cout << "exiting prematurely" << "\n";
            break;
        }
    }

    MeshTools::ReconstructMeshAfterFaceCut(m);
    MeshTools::writePLY(outputfname, m);
}

void MeshActions::ConvertToStL(CommandLineArguments& cmd){
    std::string inputfname, outputfname="out.stl";
    Mesh mesh;
    cmd.readString("i", CommandLineArguments::Keys::Required, inputfname);
    cmd.readString("o", CommandLineArguments::Keys::Optional, outputfname);
    MeshTools::readPLYlibr(inputfname, mesh);
    MeshTools::writeSTL(outputfname, mesh);
}

void MeshActions::ReplicatePeriodicMesh(CommandLineArguments& cmd){
    // box --> must provide
    std::string inputfname, outputfname="output.ply";
    Real3 box;
    INT2 translate_index={{1,2}}, array_shape;

    cmd.readString("i", CommandLineArguments::Keys::Required, inputfname);
    cmd.readString("o", CommandLineArguments::Keys::Optional, outputfname);
    cmd.readArray("box", CommandLineArguments::Keys::Required, box);
    cmd.readArray("array_shape", CommandLineArguments::Keys::Required, array_shape);
    cmd.readArray("replicate_dir_index", CommandLineArguments::Keys::Optional, translate_index);

    // read mesh
    Mesh mesh;
    MeshTools::readPLYlibr(inputfname, mesh);
    mesh.setBoxLength(box);

    // total number of replicates --> (2,2) --> 4
    int total_num_replicates = array_shape[0] * array_shape[1];

    // mesh properties
    auto& v = mesh.accessvertices();
    auto& f = mesh.accesstriangles();

    // check the periodic faces
    std::vector<int> PeriodicFaceIndices;
    std::vector<int> NonPeriodicFaceIndices;
    std::vector<triangle> newTriangles;

    // check which ones are periodic triangles 
    for (int i=0;i<f.size();i++){
        if (MeshTools::IsPeriodicTriangle(mesh, i)){
            PeriodicFaceIndices.push_back(i);
        }
        else{
            NonPeriodicFaceIndices.push_back(i);
        }
    }

    // requires periodic meshes --> replicates meshes in the 2 directions
    int t1= translate_index[0], t2=translate_index[1];
    int nv = v.size();
    std::cout << "nv = " << nv << "\n";
    std::vector<vertex> newVertices;
    for (int i=0;i<array_shape[0];i++){
        for (int j=0;j<array_shape[1];j++){
            for (int k=0;k<nv;k++){
                vertex newV = v[k];

                // shift to a new position
                newV.position_[t1] += j * box[t1];
                newV.position_[t2] += i * box[t2];
                newVertices.push_back(newV);
            }
        }
    }

    // do all the non periodic faces
    for (int i=0;i<array_shape[0];i++){
        for (int j=0;j<array_shape[1];j++){
            int num = i * array_shape[1] + j;
            for (int k=0;k<NonPeriodicFaceIndices.size();k++){
                int index = NonPeriodicFaceIndices[k];
                triangle newt = f[index];
                newt.triangleindices_ = newt.triangleindices_ + num * nv;
                newTriangles.push_back(newt);
            }
        }
    }

    // do all the periodic faces  --> rows 
    for (int i=0;i<array_shape[0];i++){
        for (int j=0;j<array_shape[1]-1;j++){
            // the numbering scheme
            /////////////////
            // 4    5    6 //
            // *    *    * //
            // 1    2    3 //
            // *    *    * //
            /////////////////
            INT2 this_ind = {{i,j}};
            int num = i * array_shape[1] + j;
            for (int m=0;m<PeriodicFaceIndices.size();m++){
                int faceIdx = PeriodicFaceIndices[m]; 

                // define the new triangle  --> first copy everything from the old triangle
                triangle newt = f[faceIdx];

                // obtain the old triangle indices --> these are real vertex indices 
                INT3 vertsIndex = f[faceIdx].triangleindices_;

                // obtain positions --> in t1 dimension
                Real3 pos = {{v[vertsIndex[0]].position_[t1], v[vertsIndex[1]].position_[t1], v[vertsIndex[2]].position_[t1]}};
                int maxIndex = Algorithm::argmax(pos);
                std::vector<int> OtherIndex;
                for (int k=0;k<3;k++){
                    if (k != maxIndex){OtherIndex.push_back(k);}
                }
                int OtherIndex1 = vertsIndex[OtherIndex[0]], OtherIndex2 = vertsIndex[OtherIndex[1]], max = vertsIndex[maxIndex];

                // find all the edge length in this triangle 
                Real3 shift1 = MeshTools::calculateShift(v[OtherIndex1].position_, v[max].position_, box);
                Real3 shift2 = MeshTools::calculateShift(v[OtherIndex2].position_, v[max].position_, box);
                INT2 shift_index1, shift_index2;
                int shift_num1=0, shift_num2=0;
                for (int k=0;k<2;k++){
                    if (shift1[translate_index[k]] == 0) {shift_index1[(k+1)%2] = 0;}
                    else if (shift1[translate_index[k]] > 0) { shift_index1[(k+1)%2] = 1;}
                    else { shift_index1[(k+1)%2] = -1;}

                    if (shift2[translate_index[k]] == 0) {shift_index2[(k+1)%2] = 0;}
                    else if (shift2[translate_index[k]] > 0) { shift_index2[(k+1)%2] = 1;}
                    else { shift_index2[(k+1)%2] = -1;}
                }
                INT2 temp1 = shift_index1 + this_ind;
                INT2 temp2 = shift_index2 + this_ind;

                if ((temp1[0] >= array_shape[0]) || (temp1[1] >= array_shape[1]) || (temp2[0] >= array_shape[0]) || (temp2[1] >= array_shape[1])){
                    continue;
                }

                shift_num1  = shift_index1[0] * array_shape[1] + shift_index1[1];
                shift_num2  = shift_index2[0] * array_shape[1] + shift_index2[1];

                int newIndex1 = shift_num1 + num;
                int newIndex2 = shift_num2 + num;
                if ((newIndex1 < 0) || (newIndex2 < 0)){
                    continue;
                }
                else if ((newIndex1 >= total_num_replicates) || (newIndex2 >= total_num_replicates)) {continue;}
                else{
                    newt.triangleindices_[maxIndex] += num * nv;
                    newt.triangleindices_[OtherIndex[0]] += newIndex1 * nv;
                    newt.triangleindices_[OtherIndex[1]] += newIndex2 * nv;
                    newTriangles.push_back(newt);
                }
            }
        }
    }

    // // do all the periodic faces  --> cols
    for (int i=0;i<array_shape[1];i++){
        for (int j=0;j<array_shape[0]-1;j++){
            INT2 this_ind = {{j,i}};
            int num = j * array_shape[1] + i;
            for (int m=0;m<PeriodicFaceIndices.size();m++){
                int faceIdx = PeriodicFaceIndices[m];

                // define the new triangle  --> first copy everything from the old triangle
                triangle newt = f[faceIdx];

                // obtain the old triangle indices --> these are real vertex indices 
                INT3 vertsIndex = f[faceIdx].triangleindices_;

                // obtain positions --> in t2 dimension
                Real3 pos = {{v[vertsIndex[0]].position_[t2], v[vertsIndex[1]].position_[t2], v[vertsIndex[2]].position_[t2]}};
                int maxIndex = Algorithm::argmax(pos);
                std::vector<int> OtherIndex;
                for (int k=0;k<3;k++){
                    if (k != maxIndex){OtherIndex.push_back(k);}
                }
                int OtherIndex1 = vertsIndex[OtherIndex[0]], OtherIndex2 = vertsIndex[OtherIndex[1]], max = vertsIndex[maxIndex];

                // find all the edge length in this triangle 
                Real3 shift1 = MeshTools::calculateShift(v[OtherIndex1].position_, v[max].position_, box);
                Real3 shift2 = MeshTools::calculateShift(v[OtherIndex2].position_, v[max].position_, box);
                INT2 shift_index1, shift_index2;
                int shift_num1=0, shift_num2=0;
                for (int k=0;k<2;k++){
                    if (shift1[translate_index[k]] == 0) {shift_index1[(k+1)%2] = 0;}
                    else if (shift1[translate_index[k]] > 0) { shift_index1[(k+1)%2] = 1;}
                    else { shift_index1[(k+1)%2] = -1;}

                    if (shift2[translate_index[k]] == 0) {shift_index2[(k+1)%2] = 0;}
                    else if (shift2[translate_index[k]] > 0) { shift_index2[(k+1)%2] = 1;}
                    else { shift_index2[(k+1)%2] = -1;}
                }

                INT2 temp1 = shift_index1 + this_ind;
                INT2 temp2 = shift_index2 + this_ind;

                if ((temp1[0] >= array_shape[0]) || (temp1[1] >= array_shape[1]) || (temp2[0] >= array_shape[0]) || (temp2[1] >= array_shape[1])){
                    continue;
                }
                
                shift_num1  = shift_index1[0] * array_shape[1] + shift_index1[1];
                shift_num2  = shift_index2[0] * array_shape[1] + shift_index2[1];

                int newIndex1 = shift_num1 + num;
                int newIndex2 = shift_num2 + num;
                if ((newIndex1 < 0) || (newIndex2 < 0)){
                    continue;
                }
                else if ((newIndex1 >= total_num_replicates) || (newIndex2 >= total_num_replicates)) {continue;}
                else{
                    newt.triangleindices_[maxIndex] += num * nv;
                    newt.triangleindices_[OtherIndex[0]] += newIndex1 * nv;
                    newt.triangleindices_[OtherIndex[1]] += newIndex2 * nv;
                    
                    newTriangles.push_back(newt);
                }
            }
        }
    }

    v.clear();
    v.insert(v.end(), newVertices.begin(), newVertices.end());
    f.clear();
    f.insert(f.end(), newTriangles.begin(), newTriangles.end());

    std::vector<Real3> V;
    std::vector<INT3> F;
    MeshTools::RemoveDuplicatedFaces(mesh);
    MeshTools::RemoveIsolatedFaces(mesh);
    MeshTools::ConvertToNonPBCMesh(mesh);
    MeshTools::RemoveIsolatedVertices(mesh);
    mesh.CalcVertexNormals();
    MeshTools::writePLY(outputfname, mesh);
}

void MeshActions::TriangleAngleDistribution(CommandLineArguments& cmd){
    using Range = CommonTypes::Real2;
    using Binptr= std::unique_ptr<Bin>;

    std::string inputfname, outputfname="dist.out";
    Real3 box;
    Real angle=120;
    int numbins;

    cmd.readValue("i", CommandLineArguments::Keys::Required, inputfname);
    cmd.readValue("o", CommandLineArguments::Keys::Optional, outputfname);
    cmd.readValue("numbins", CommandLineArguments::Keys::Optional, numbins);
    bool BoxRead = cmd.readArray("box", CommandLineArguments::Keys::Optional, box);

    Mesh m;
    MeshTools::readPLYlibr(inputfname, m);
    if (BoxRead){
        m.setBoxLength(box);
    }

    // create the bin object
    Range r = {{0,180}};
    Bin b(r, numbins);

    // specify all the angles 
    std::vector<Real> angles;
    MeshTools::FindTriangleAngles(m, angles);

    // define the histogram
    std::vector<int> histogram_angles(b.getNumbins(), 0);

    #pragma omp parallel
    {
        std::vector<int> local_histogram(b.getNumbins(),0);
        #pragma omp for
        for (int i=0;i<angles.size();i++){
            if (b.isInRange(angles[i])){
                int binNum = b.findBin(angles[i]);
                local_histogram[binNum] += 1;
            }
        }

        #pragma omp critical
        {
            histogram_angles = histogram_angles + local_histogram;
        }
    }

    std::ofstream ofs;
    ofs.open(outputfname);

    for (int i=0;i<histogram_angles.size();i++){
        ofs << histogram_angles[i] << "\n";
    }

    ofs.close();
}

void MeshActions::SideLengthsDistribution(CommandLineArguments& cmd){
    std::string inputfname, outputfname="SideLengthDist.out";
    Real3 box;

    cmd.readValue("i", CommandLineArguments::Keys::Required, inputfname);
    cmd.readValue("o", CommandLineArguments::Keys::Optional, outputfname);
    bool ispbc = cmd.readArray("box", CommandLineArguments::Keys::Optional, box);

    Mesh m;
    MeshTools::readPLYlibr(inputfname, m);
    if (ispbc){
        m.setBoxLength(box);
    }

    std::vector<Real> SideLengths;
    MeshTools::FindSideLengths(m, SideLengths);

    std::ofstream ofs;
    ofs.open(outputfname);
    for (int i=0;i<SideLengths.size();i++){
        ofs << SideLengths[i] << "\n";
    }
    ofs.close();
}

void MeshActions::MeshCleanup(CommandLineArguments& cmd){
    using INT2 = CommonTypes::index2;
    std::string inputfname, outputfname="cleaned.ply", fixed_index_file;

    Real3 box;
    Real edgeLength, angleThreshold=120;
    int maxiterations=10, min_numneighbors=0, search_ring=5;
    bool verbose=false, preserve_features=true;
    cmd.readValue("i", CommandLineArguments::Keys::Required, inputfname);
    cmd.readValue("maxiteration", CommandLineArguments::Keys::Optional, maxiterations);
    cmd.readValue("o", CommandLineArguments::Keys::Optional, outputfname);
    cmd.readValue("angleThreshold", CommandLineArguments::Keys::Optional, angleThreshold);
    cmd.readBool("verbose", CommandLineArguments::Keys::Optional, verbose);
    cmd.readBool("preserve_features", CommandLineArguments::Keys::Optional, preserve_features);
    cmd.readValue("search_ring", CommandLineArguments::Keys::Optional, search_ring);
    cmd.readValue("min_numneighbors", CommandLineArguments::Keys::Optional, min_numneighbors);
    bool fixed = cmd.readString("fixed_index_file", CommandLineArguments::Keys::Optional, fixed_index_file);
    bool edgeLengthRead = cmd.readValue("edgeLengthCutoff", CommandLineArguments::Keys::Optional, edgeLength);
    bool isPBC = cmd.readArray("box", CommandLineArguments::Keys::Optional, box);

    Mesh m;
    MeshTools::readPLYlibr(inputfname, m);
    if (isPBC){
        m.setBoxLength(box);
    }

    int nv = m.getNumVertices();

    std::vector<int> fixed_index;
    std::vector<int> isfixed(nv, 0);
    if (fixed){
        StringTools::ReadTabulatedData(fixed_index_file, 0, fixed_index);
        for (int i=0;i<fixed_index.size();i++){
            isfixed[fixed_index[i]] = 1;
        }
    }

    if (min_numneighbors > 0){
        MeshTools::RemoveMinimumNeighbors(m, search_ring, min_numneighbors);
    }

    // if edge length is not provided, then we use the average side length
    if (! edgeLengthRead){
        std::vector<Real> sideLength;
        MeshTools::FindSideLengths(m, sideLength);
        Real sum=0.0;
        #pragma omp parallel
        {
            Real localsum = 0.0;
            #pragma omp for
            for (int i=0;i<sideLength.size();i++){
                localsum  += sideLength[i];
            }

            #pragma omp critical
            {
                sum += localsum;
            }
        }

        edgeLength = sum / sideLength.size();

        if (verbose){
            std::cout << "edge length cut off  = " << edgeLength << "\n";
        }
    }

    int count=0;
    // initialize edge removal 



    while (true){
        int verticesbefore = m.getvertices().size();

        // first remove short edges
        ShortEdgeRemoval edgeRemove(m);

        // if we are preserving features, then we set high importance for the boundary vertices
        std::vector<int> importance(verticesbefore,0);
        if (preserve_features){
            std::map<INT2,std::vector<int>> EToF; 
            std::vector<bool> boundaryIndicator;
            std::vector<int> importance(verticesbefore,0);
            MeshTools::MapEdgeToFace(m, EToF);
            MeshTools::CalculateBoundaryVertices(m, EToF, boundaryIndicator);
            for (int i=0;i<verticesbefore;i++){
                if (MeshTools::IsBoundary(i, boundaryIndicator)){
                    // high importance 
                    importance[i] = -1;
                }
            }
        }

        if (fixed){
            for (int i=0;i<verticesbefore;i++){
                if (isfixed[i]){
                    importance[i] = -1;
                }
            }
        }

        edgeRemove.set_importance(importance);

        int numCollapsed = edgeRemove.calculate(edgeLength);

        // then remove obtuse edges 
        ObtuseTriangleRemoval triangleRemove(m);
        triangleRemove.calculate(angleThreshold,100);

        int facesafter = m.gettriangles().size(); 
        if (verbose){
            std::cout << "Collapsed edges = " << numCollapsed << " facesAfter = " << facesafter << "\n";
        }
        int verticesafter = m.getvertices().size();
        if (verticesafter == verticesbefore){
            break;
        }

        if (count > maxiterations){
            std::cout << "The clean up process did not converge, but exit before it exceed max iteration number " << maxiterations << "\n";
            break;
        }
        MeshTools::RemoveIsolatedVertices(m);
        MeshTools::RemoveDuplicatedFaces(m);
        MeshTools::RemoveIsolatedFaces(m);

        count ++;
    }

    MeshTools::RemoveIsolatedVertices(m);
    MeshTools::RemoveDuplicatedFaces(m);
    MeshTools::RemoveIsolatedFaces(m);


    m.CalcVertexNormals();


    MeshTools::writePLY(outputfname, m);
}

void MeshActions::calculateSurfaceProperties(CommandLineArguments& cmd){
    // declare the shape 
    std::unique_ptr<AFP_shape> shape;
    std::string inputfname, outputfname;
    std::string surface_shape_name;
    ParameterPack param;

    cmd.readString("i", CommandLineArguments::Keys::Required, inputfname);
    cmd.readString("o", CommandLineArguments::Keys::Optional, outputfname);
    cmd.readString("shape", CommandLineArguments::Keys::Required, surface_shape_name);

    if (surface_shape_name == "superegg"){
        std::string a,b,n,zmax,a_taper,b_taper,a_alpha,b_alpha;
        std::vector<std::string> center;
        bool read;
        cmd.readValue("a", CommandLineArguments::Keys::Required, a);
        param.insert("a", a);
        cmd.readValue("b", CommandLineArguments::Keys::Required, b);
        param.insert("b", b);
        cmd.readValue("zmax", CommandLineArguments::Keys::Required, zmax);
        param.insert("zmax", zmax);
        cmd.readVector("center", CommandLineArguments::Keys::Required, center);
        param.insert("center", center);
        read = cmd.readValue("n", CommandLineArguments::Keys::Optional, n);
        if (read){
            param.insert("n", n);
        }
        read = cmd.readValue("a_taper", CommandLineArguments::Keys::Optional, a_taper);
        if (read){
            param.insert("a_taper", a_taper);
        }
        read = cmd.readValue("b_taper", CommandLineArguments::Keys::Optional, b_taper);
        if (read){
            param.insert("b_taper", b_taper);
        }
        read = cmd.readValue("a_alpha", CommandLineArguments::Keys::Optional, a_alpha);
        if (read){
            param.insert("a_alpha", a_alpha);
        }
        read = cmd.readValue("b_alpha", CommandLineArguments::Keys::Optional, b_alpha);
        if (read){
            param.insert("b_alpha", b_alpha);
        }

        // initialize the shape 
        shape = std::make_unique<SuperEgg>(param);
    }
    else{
        ASSERT((true == false), "The shape " << surface_shape_name << " is not implemented yet.");
    }

    Real num_Steps = 100;
    Real v_step = Constants::PI / (2 * num_Steps);
    std::vector<Real3> normals, tangents;
    std::vector<double3> pos;

    for (float v = 0; v < Constants::PI / 2; v += v_step){
        Real3 tangent, normal;
        shape->CalculateNormalAndTangent(0, v, tangent, normal, 0, 1, 2);
        pos.push_back(shape->calculatePos(0,v));
        normals.push_back(normal);
        tangents.push_back(tangent);
    }

    std::ofstream ofs;
    ofs.open("output.out");
    for (int i=0;i<normals.size();i++){
        ofs << pos[i][0] << " " << pos[i][1] << " " << pos[i][2] << " " <<  \
               normals[i][0] << " " << normals[i][1] << " " << normals[i][2] << \
         " " << tangents[i][0] << " " << tangents[i][1] << " " << tangents[i][2] << "\n";
    }
    ofs.close();
}

    
void MeshActions::FlattenMesh(CommandLineArguments& cmd){
    std::string inputfname, outputfname="flat.ply";
    int index=2;

    cmd.readString("i", CommandLineArguments::Keys::Required, inputfname);
    cmd.readString("o", CommandLineArguments::Keys::Optional, outputfname);
    cmd.readValue("index", CommandLineArguments::Keys::Required, index);

    Mesh m;
    MeshTools::readPLYlibr(inputfname, m);

    auto& vertices = m.accessvertices();
    for (int i=0;i<vertices.size();i++){
        vertices[i].position_[index] = 0;
    }

    MeshTools::writePLY(outputfname, m);
}

void MeshActions::ConformingTriangulations(CommandLineArguments& cmd){
    std::string inputfname, outputfname="gen.ply", boundaryfname;
    std::vector<std::string> boundaryfnames;
    Real2 Box;
    Real aspect_bound=0.125, size_bound=0.5;
    int NumBoundary;
    Real height=1.0;
    bool periodic=true;
    bool isBoundingBox=true;
    INT2 index = {{0,1}};

    cmd.readString("o", CommandLineArguments::Keys::Optional, outputfname);
    cmd.readVector("boundaryfiles", CommandLineArguments::Keys::Required, boundaryfnames);
    cmd.readBool("isBoundingBox", CommandLineArguments::Keys::Optional, isBoundingBox);
    cmd.readValue("aspect_bound", CommandLineArguments::Keys::Optional, aspect_bound);
    cmd.readValue("size_bound", CommandLineArguments::Keys::Optional, size_bound);
    cmd.readBool("periodic", CommandLineArguments::Keys::Optional, periodic);
    cmd.readValue("height", CommandLineArguments::Keys::Optional, height);

    if (periodic){
        cmd.readArray("box", CommandLineArguments::Keys::Required, Box);
    }
    cmd.readArray("ReadIndex", CommandLineArguments::Keys::Optional, index);

    if (isBoundingBox){
        cmd.readValue("NumBoundary", CommandLineArguments::Keys::Required, NumBoundary);
    }

    // copy the boundary points 
    std::vector<Real2> points;
    std::vector<INT2> edges;
    std::vector<Real> edgeLength;
    std::vector<Real2> seed;

    for (auto f : boundaryfnames){
        // define a variable that holds the centroid of the shape
        Real2 centroid={{0,0}};

        // a temporary vector to store data
        std::vector<std::vector<Real>> temp;
        StringTools::ReadTabulatedData(f, temp);
        ASSERT((temp.size() != 0), "The boundary file provided contains no data.");

        // write the data in temp into points 
        int points_initial_size = points.size();
        for (int i=0;i<temp.size();i++){
            Real2 pos;
            INT2 e;
            pos[0] = temp[i][index[0]];
            pos[1] = temp[i][index[1]];
            centroid = centroid + pos;
            points.push_back(pos);
        }

        // find the centroid of the points 
        centroid = centroid / temp.size();
        if (isBoundingBox){
            seed.push_back(centroid);
        }

        // get rid of bad edges
        for (int i=points_initial_size;i<points.size();i++){
            int next = (i + 1) % points.size();
            if (next == 0){
                next += points_initial_size;
            }
            Real2 diff = points[i] - points[next];
            Real sum=0;
            for (int j=0;j<2;j++){
                sum += diff[j] * diff[j];
            }
            if (sum > 1e-6){
                INT2 e = {{i,next}};
                edges.push_back(e);
                edgeLength.push_back(std::sqrt(sum));
            }
        }

        Real minEdgeLength = Algorithm::min(edgeLength);
    }
    std::cout << edges << std::endl;

    if (isBoundingBox){
        // construct the box
        Real2 stepsize = Box / NumBoundary;
        std::vector<Real2> newPoints;
        std::vector<INT2> newEdges;
        newPoints.push_back({{0,0}});
        // start with x = 0 y = something
        for (int i=1;i<NumBoundary+1;i++){
            Real2 p = {{0,i * stepsize[1]}};
            newPoints.push_back(p);
            newEdges.push_back({{i-1,i}});
        }

        // continue with x = something, y = max
        for (int i=1;i<NumBoundary+1;i++){
            Real2 p = {{i * stepsize[0], Box[1]}};
            int LastIndex = newPoints.size() - 1;
            newPoints.push_back(p);
            int CurrIndex = newPoints.size() - 1;
            newEdges.push_back({{LastIndex, CurrIndex}});
        }

        // continnue with x = max ,  y =something
        for (int i=NumBoundary-1;i>=0;i--){
            Real2 p = {{Box[0], i * stepsize[1]}};
            int LastIndex = newPoints.size()-1;
            newPoints.push_back(p);
            int CurrIndex = newPoints.size() - 1;
            newEdges.push_back({{LastIndex, CurrIndex}});
        }

        // continue with x = something, y = 0
        for (int i=NumBoundary-1;i>=1;i--){
            Real2 p = {{i * stepsize[0], 0}};
            int LastIndex = newPoints.size() - 1;
            newPoints.push_back(p);
            int CurrIndex = newPoints.size() - 1;
            newEdges.push_back({{LastIndex, CurrIndex}});
        }

        // finally connect the last point with the first point
        int LastIndex = newPoints.size()-1;
        newEdges.push_back({{LastIndex, 0}});
        int ori_point_size = points.size();
        newEdges = newEdges + ori_point_size;

        // merge the newpoints with points 
        points.insert(points.end(), newPoints.begin(), newPoints.end());
        edges.insert(edges.end(), newEdges.begin(), newEdges.end());
    }

    // constrct the generation 
    MeshGen2d meshgen(points, edges, seed, aspect_bound, size_bound, periodic);
    if (periodic){
        meshgen.setBoxLength(Box);
    }
    meshgen.generate();
    const auto& m = meshgen.getMesh();
    MeshTools::writePLY(outputfname, m);
}

void MeshActions::SplitLongEdges(CommandLineArguments& cmd){
    std::string inputfname, outputfname="split.ply";

    bool isPBC;
    Real3 box;
    Real max_length;
    cmd.readString("i", CommandLineArguments::Keys::Required, inputfname);
    cmd.readString("o", CommandLineArguments::Keys::Optional, outputfname);
    isPBC = cmd.readArray("box", CommandLineArguments::Keys::Optional, box);
    cmd.readValue("max_length", CommandLineArguments::Keys::Required, max_length);

    Mesh m;
    MeshTools::readPLYlibr(inputfname, m);

    if (isPBC){
        m.setBoxLength(box);
    }

    LongEdgeRemoval e(m);
    e.run(max_length);

    MeshTools::RemoveIsolatedVertices(m);
    MeshTools::RemoveDuplicatedFaces(m);
    MeshTools::RemoveIsolatedFaces(m);
    m.CalcVertexNormals();

    MeshTools::writePLY(outputfname, m);
}

void MeshActions::ShiftMeshWithRef(CommandLineArguments& cmd){
    std::string inputfname, outputfname="shifted.ply", reffname;

    cmd.readString("i", CommandLineArguments::Keys::Required, inputfname);
    cmd.readString("o", CommandLineArguments::Keys::Optional, outputfname);
    cmd.readString("ref", CommandLineArguments::Keys::Required, reffname);

    Mesh m, ref;
    MeshTools::readPLYlibr(inputfname, m);
    MeshTools::readPLYlibr(reffname, ref);

    auto& v = m.accessvertices();
    const auto& refv = ref.getvertices();

    Real3 v_COM = {{0,0,0}};
    Real3 ref_COM = {{0,0,0}};
    for (int i=0;i<v.size();i++){
        v_COM = v_COM + v[i].position_;
    }

    for (int i=0;i<refv.size();i++){
        ref_COM = ref_COM + refv[i].position_;
    }

    v_COM = v_COM / v.size();
    ref_COM = ref_COM / refv.size();
    Real3 diff = ref_COM - v_COM;

    for (int i=0;i<v.size();i++){
        v[i].position_ = v[i].position_ + diff;
    }


    MeshTools::writePLY(outputfname, m);
}


void MeshActions::ViewMeshWithData(CommandLineArguments& cmd){
    #ifdef IGL_ENABLED
    std::string inputfname, datafname;
    int col, numSteps=21;
    double min, max;
    Real3 box;
    cmd.readValue("i", CommandLineArguments::Keys::Required, inputfname);
    cmd.readValue("data", CommandLineArguments::Keys::Required, datafname);
    cmd.readValue("col", CommandLineArguments::Keys::Required, col);
    cmd.readValue("numSteps", CommandLineArguments::Keys::Optional, numSteps);
    bool minRead = cmd.readValue("min", CommandLineArguments::Keys::Optional,min); 
    bool maxRead = cmd.readValue("max", CommandLineArguments::Keys::Optional, max);
    bool pbc = cmd.readArray("box", CommandLineArguments::Keys::Optional, box);

    // mesh read with my program
    Mesh m;
    MeshTools::readPLYlibr(inputfname,m);
    if (pbc){
        m.setBoxLength(box);
    }

    // mesh read with igl
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    igl::read_triangle_mesh(inputfname, V, F);

    std::vector<double> data;
    StringTools::ReadTabulatedData(datafname, col, data);

    if (pbc){
        const auto& tri = m.gettriangles();
        const auto& v = m.getvertices();
        std::vector<vertex> newv;
        newv.insert(newv.end(), v.begin(), v.end());

        std::vector<int> nonPeriodicFace;
        std::vector<triangle> newT;
        for (int i=0;i<tri.size();i++){
            if (MeshTools::IsPeriodicTriangle(m,i)){
                // shift with respect to A
                Real3 AA,BA,CA, AB,BB,CB, AC,BC,CC, centroid_refA, centroid_refB, centroid_refC;
                Real3 Aref, Bref, Cref;
                Real d1,d2,d3;
                Aref = v[tri[i][0]].position_;
                Bref = v[tri[i][1]].position_;
                Cref = v[tri[i][2]].position_;
                MeshTools::ShiftPeriodicTriangle(v, tri[i].triangleindices_, m.getBoxLength(),Aref, AA,BA,CA);
                MeshTools::ShiftPeriodicTriangle(v, tri[i].triangleindices_, m.getBoxLength(),Bref, AB,BB,CB);
                MeshTools::ShiftPeriodicTriangle(v, tri[i].triangleindices_, m.getBoxLength(),Cref, AC,BC,CC);
                d1 = data[tri[i][0]]; d2 = data[tri[i][1]]; d3 = data[tri[i][2]];

                centroid_refA = (AA+BA+CA) / 3.0;
                centroid_refB = (AB+BB+CB) / 3.0;
                centroid_refC = (AC+CB+CC) / 3.0;

                // append A
                int index = newv.size();                
                vertex v1,v2,v3;
                v1.position_ = AA; v2.position_=BA; v3.position_=CA;
                newv.push_back(v1); newv.push_back(v2); newv.push_back(v3);
                INT3 tIndex = {{index,index+1,index+2}};
                triangle t;
                t.triangleindices_ = tIndex;
                newT.push_back(t);
                data.push_back(d1); data.push_back(d2); data.push_back(d3);

                // append B
                Real3 diffAB = centroid_refB - centroid_refA;
                Real diffsq_AB=0.0;
                for (int j=0;j<3;j++){
                    diffsq_AB += (diffAB[i] * diffAB[i]);
                }
                diffsq_AB = std::sqrt(diffsq_AB);
                if (diffsq_AB < 1e-5){
                    int index = newv.size();                
                    v1.position_ = AB; v2.position_=BB; v3.position_=CB;
                    newv.push_back(v1); newv.push_back(v2); newv.push_back(v3);
                    INT3 tIndex = {{index,index+1,index+2}};
                    t.triangleindices_ = tIndex;
                    newT.push_back(t);
                    data.push_back(d1); data.push_back(d2); data.push_back(d3);
                }

                // append C
                Real3 diffCB = centroid_refC - centroid_refB;
                Real diffsq_CB=0.0;
                for (int j=0;j<3;j++){
                    diffsq_CB += (diffCB[i] * diffCB[i]);
                }
                diffsq_CB = std::sqrt(diffsq_CB);
                if (diffsq_CB < 1e-5){
                    int index = newv.size();                
                    v1.position_ = AC; v2.position_=BC; v3.position_=CC;
                    newv.push_back(v1); newv.push_back(v2); newv.push_back(v3);
                    INT3 tIndex = {{index,index+1,index+2}};
                    triangle t;
                    t.triangleindices_ = tIndex;
                    newT.push_back(t);
                    data.push_back(d1); data.push_back(d2); data.push_back(d3);
                }
            }
            else{
                newT.push_back(tri[i]);
            }

        }

        Eigen::MatrixXi temp(newT.size(),3);
        for (int i=0;i<newT.size();i++){
            for (int j=0;j<3;j++){
                temp(i,j) = newT[i].triangleindices_[j];
            }
        }

        Eigen::MatrixXd newM(newv.size(),3);
        for (int i=0;i<newv.size();i++){
            for (int j=0;j<3;j++){
                newM(i,j) = newv[i].position_[j];

            }
        }

        V = newM;
        F = temp;
    }
    
    Eigen::VectorXd d = Eigen::Map<Eigen::VectorXd>(data.data(), data.size());
    ASSERT((d.rows() == V.rows()), "The data size " << d.rows() << " does not match the vertics size " << V.rows());

    if (! minRead){min = Algorithm::min(data);}
    if (! maxRead){max = Algorithm::max(data);}

    // set up viewer
    igl::opengl::glfw::Viewer viewer;
    viewer.data().set_mesh(V,F);
    viewer.data().set_data(d.transpose(), min, max, igl::ColorMapType::COLOR_MAP_TYPE_JET, numSteps);
    viewer.data().show_lines = false;
    viewer.launch();
    #else
    std::cout << "IGL is not enabled." << "\n";
    #endif
}

void MeshActions::ViewMesh(CommandLineArguments& cmd){
    #ifdef IGL_ENABLED
    std::string inputfname;
    Real3 box;
    cmd.readValue("i", CommandLineArguments::Keys::Required, inputfname);
    bool pbc = cmd.readArray("box", CommandLineArguments::Keys::Optional, box);

    // mesh read with my program
    Mesh m;
    MeshTools::readPLYlibr(inputfname,m);
    if (pbc){
        m.setBoxLength(box);
    }

    // mesh read with igl
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    igl::readPLY(inputfname, V, F);

    Eigen::MatrixXd N_vertices, N_faces, N_corners; 

    // This function is called every time a keyboard button is pressed
    auto key_down = [N_vertices, N_faces, N_corners](igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifier) -> bool
    {
        switch(key)
        {
            case '1':
            viewer.data().set_normals(N_faces);
            return true;
            case '2':
            viewer.data().set_normals(N_vertices);
            return true;
            case '3':
            viewer.data().set_normals(N_corners);
            return true;
            default: break;
        }
    return false;
    };

    if (pbc){
        const auto& tri = m.gettriangles();
        const auto& v = m.getvertices();
        std::vector<vertex> newv;
        newv.insert(newv.end(), v.begin(), v.end());

        std::vector<int> nonPeriodicFace;
        std::vector<triangle> newT;
        for (int i=0;i<tri.size();i++){
            if (MeshTools::IsPeriodicTriangle(m,i)){
                // shift with respect to A
                Real3 AA,BA,CA, AB,BB,CB, AC,BC,CC, centroid_refA, centroid_refB, centroid_refC;
                Real3 Aref, Bref, Cref;
                Real d1,d2,d3;
                Aref = v[tri[i][0]].position_;
                Bref = v[tri[i][1]].position_;
                Cref = v[tri[i][2]].position_;
                MeshTools::ShiftPeriodicTriangle(v, tri[i].triangleindices_, m.getBoxLength(),Aref, AA,BA,CA);
                MeshTools::ShiftPeriodicTriangle(v, tri[i].triangleindices_, m.getBoxLength(),Bref, AB,BB,CB);
                MeshTools::ShiftPeriodicTriangle(v, tri[i].triangleindices_, m.getBoxLength(),Cref, AC,BC,CC);

                centroid_refA = (AA+BA+CA) / 3.0;
                centroid_refB = (AB+BB+CB) / 3.0;
                centroid_refC = (AC+CB+CC) / 3.0;

                // append A
                int index = newv.size();                
                vertex v1,v2,v3;
                v1.position_ = AA; v2.position_=BA; v3.position_=CA;
                newv.push_back(v1); newv.push_back(v2); newv.push_back(v3);
                INT3 tIndex = {{index,index+1,index+2}};
                triangle t;
                t.triangleindices_ = tIndex;
                newT.push_back(t);

                // append B
                Real3 diffAB = centroid_refB - centroid_refA;
                Real diffsq_AB=0.0;
                for (int j=0;j<3;j++){
                    diffsq_AB += (diffAB[i] * diffAB[i]);
                }
                diffsq_AB = std::sqrt(diffsq_AB);
                if (diffsq_AB < 1e-5){
                    int index = newv.size();                
                    v1.position_ = AB; v2.position_=BB; v3.position_=CB;
                    newv.push_back(v1); newv.push_back(v2); newv.push_back(v3);
                    INT3 tIndex = {{index,index+1,index+2}};
                    t.triangleindices_ = tIndex;
                    newT.push_back(t);
                }

                // append C
                Real3 diffCB = centroid_refC - centroid_refB;
                Real diffsq_CB=0.0;
                for (int j=0;j<3;j++){
                    diffsq_CB += (diffCB[i] * diffCB[i]);
                }
                diffsq_CB = std::sqrt(diffsq_CB);
                if (diffsq_CB < 1e-5){
                    int index = newv.size();                
                    v1.position_ = AC; v2.position_=BC; v3.position_=CC;
                    newv.push_back(v1); newv.push_back(v2); newv.push_back(v3);
                    INT3 tIndex = {{index,index+1,index+2}};
                    triangle t;
                    t.triangleindices_ = tIndex;
                    newT.push_back(t);
                }
            }
            else{
                newT.push_back(tri[i]);
            }

        }

        Eigen::MatrixXi temp(newT.size(),3);
        for (int i=0;i<newT.size();i++){
            for (int j=0;j<3;j++){
                temp(i,j) = newT[i].triangleindices_[j];
            }
        }

        Eigen::MatrixXd newM(newv.size(),3);
        for (int i=0;i<newv.size();i++){
            for (int j=0;j<3;j++){
                newM(i,j) = newv[i].position_[j];

            }
        }

        V = newM;
        F = temp;
    }


    igl::opengl::glfw::Viewer viewer;
    igl::per_face_normals(V,F,N_faces);
    igl::per_vertex_normals(V,F,N_vertices);
    viewer.callback_key_down = key_down;
      std::cout<<
    "Press '1' for per-face normals."<<std::endl<<
    "Press '2' for per-vertex normals."<<std::endl<<
    "Press '3' for per-corner normals."<<std::endl;
    viewer.data().set_mesh(V,F);
    viewer.data().set_normals(N_faces);
    viewer.data().show_lines = false;

    int err = viewer.launch();
    #else
    std::cout << "IGL is not enabled." << "\n";
    #endif
}



void MeshActions::FlattenMeshDimension(CommandLineArguments& cmd)
{
    std::string inputfname, outputfname="flat.ply";
    int dim;
    Real set_value=0.0;
    cmd.readString("i", CommandLineArguments::Keys::Required, inputfname);
    cmd.readString("o", CommandLineArguments::Keys::Optional, outputfname);
    cmd.readValue("dim", CommandLineArguments::Keys::Required, dim);
    cmd.readValue("set", CommandLineArguments::Keys::Optional, set_value);

    Mesh m;
    MeshTools::readPLYlibr(inputfname, m);

    auto& v = m.accessvertices();

    for (int i=0;i<v.size();i++){
        v[i].position_[dim] = set_value;
    }

    MeshTools::writePLY(outputfname, m);
}

void MeshActions::DistanceBetweenMeshesMT(CommandLineArguments& cmd){
    std::array<std::string,2> inputfnames;
    std::string outputfname;
    Real3 box, RayDirection, oppositeRay;

    // read the inputs
    cmd.readArray("i", CommandLineArguments::Keys::Required,inputfnames);
    cmd.readString("o", CommandLineArguments::Keys::Optional, outputfname);
    cmd.readArray("RayDirection", CommandLineArguments::Keys::Required, RayDirection);
    oppositeRay = -1 * RayDirection;
    bool isPBC = cmd.readArray("box", CommandLineArguments::Keys::Optional, box);


    Mesh m1, m2;
    MeshTools::readPLYlibr(inputfnames[0], m1); MeshTools::readPLYlibr(inputfnames[1], m2);

    // set the box length if is pbc
    if (isPBC){
        m1.setBoxLength(box);
        m2.setBoxLength(box);
    }

    // obtain the vertices and triangles 
    const auto& v1 = m1.getvertices();
    const auto& f1 = m1.gettriangles();
    const auto& v2 = m2.getvertices();
    const auto& f2 = m2.gettriangles();

    // first we check the boundary vertices --> if boundary vertices, then we don't calculate the difference 
    std::map<INT2, std::vector<int>> MapEdgeToFace;
    std::vector<bool> boundaryIndicator;
    MeshTools::MapEdgeToFace(m1, MapEdgeToFace);
    MeshTools::CalculateBoundaryVertices(m1, MapEdgeToFace, boundaryIndicator);

    std::vector<Real> VecDist(v1.size(),-1);

    #pragma omp parallel for
    for (int i=0;i<v1.size();i++){
        Real minDist = std::numeric_limits<Real>::max();
        Real3 O = v1[i].position_;
        bool hitOnce = false;
        if (! boundaryIndicator[i]){
            for (int j=0;j<f2.size();j++){
                Real3 A,B,C;
                if (isPBC){
                    MeshTools::ShiftPeriodicTriangle(v2, f2[j].triangleindices_, box, v1[i].position_, A, B, C);
                }
                else{
                    A = v2[f2[j][0]].position_; B = v2[f2[j][1]].position_; C = v2[f2[j][2]].position_;
                }

                Real t,u,v;

                // see if Ray Direction hits 
                bool hit = MeshTools::MTRayTriangleIntersection(A,B,C, O, RayDirection, t, u, v);
                // if not try the opposite ray
                if (! hit){
                    hit = MeshTools::MTRayTriangleIntersection(A,B,C,O,oppositeRay,t,u,v);
                }

                if (hit && std::abs(t) < minDist){
                    minDist = std::abs(t);
                    hitOnce = true;
                }
            }

            ASSERT((hitOnce), "The vertex " << O << " did not hit any triangles.");
            VecDist[i] = minDist;
        }
    }

    std::ofstream ofs;
    ofs.open(outputfname);

    for (Real vec : VecDist){
        ofs << vec << "\n";
    }

    ofs.close();
}

void MeshActions::IterativeClosestPoint(CommandLineArguments& cmd){
    std::array<std::string,2> inputfnames;
    Real3 box;

    cmd.readArray("i", CommandLineArguments::Keys::Required, inputfnames);
    bool isPBC = cmd.readArray("box", CommandLineArguments::Keys::Optional, box);

    Mesh m1, m2;
    MeshTools::readPLYlibr(inputfnames[0], m1); MeshTools::readPLYlibr(inputfnames[1], m2);

    // if (isPBC){
    //     m1.setBoxLength(box);
    //     m2.setBoxLength(box);
    // }

    const auto& v1 = m1.getvertices();
    const auto& v2 = m2.getvertices();

    //ASSERT((v1.size() == v2.size()), "The number of vertex must be equal for now.");
    std::vector<Real3> pos1(v1.size()), pos2(v2.size());

    for (int i=0;i<v1.size();i++){
        pos1[i] = v1[i].position_; pos2[i] = v2[i].position_;
    }

    Real3 COM1={0,0,0}, COM2={0,0,0};
    for (int i=0;i<pos2.size();i++){
        COM2 = COM2 + pos2[i];
    }

    COM2 = COM2 / pos2.size();

    for (int i=0;i<pos2.size();i++){
        pos2[i] = pos2[i] - COM2;
    }

    int iteration=0;
    while (true){
        // first subtract the COM
        COM1={0,0,0};
        for (int i=0;i<pos1.size();i++){
            COM1 = COM1 + pos1[i];
        }

        COM1 = COM1 / pos1.size();

        // subtract COM1 from pos1, COM2 from pos2
        std::vector<Real3> pos1_centered(pos1.size());
        for (int i=0;i<pos1.size();i++){
            pos1_centered[i] = pos1[i] - COM1;
        }

        // then find some kind of correspondence between pos1 and pos2
        std::vector<int> Association(pos1.size(),-1);
        Real total_dist=0.0;
        #pragma omp parallel
        {
            Real local_dist=0.0;
            #pragma omp for
            for (int i=0;i<pos1_centered.size();i++){
                Real min = std::numeric_limits<Real>::max();
                int minIndex = -1;
                for (int j=0;j<pos2.size();j++){
                    Real3 distVec;
                    Real dist;
                    m1.getVertexDistance(pos1_centered[i], pos2[j], distVec, dist);

                    if (dist < min){
                        minIndex = j;
                        min = dist;
                    }
                }

                ASSERT((minIndex != -1), "Something went wrong.");
                Association[i] = minIndex;
                local_dist += min;
            }

            #pragma omp critical
            {
                total_dist += local_dist;
            }
        }

        Real RMSD = total_dist / pos1.size();
        std::cout << "RMSD = " << RMSD << "\n";

        // then calculate covariance
        LinAlg3x3::Matrix cov = {};
        for (int i=0;i<Association.size();i++){
            int associated_index = Association[i];
            auto mat = LinAlg3x3::dyad(pos1_centered[i], pos2[associated_index]);

            cov = cov + mat;
        }

        Eigen::Matrix3f m;
        m << cov[0][0], cov[0][1], cov[0][2], cov[1][0], cov[1][1], cov[1][2], cov[2][0], cov[2][1], cov[2][2];
        Eigen::JacobiSVD<Eigen::Matrix3f> svd(m, Eigen::ComputeFullU | Eigen::ComputeFullV);
        auto U = svd.matrixU();
        auto V = svd.matrixV();

        auto R_found = U * V.transpose();
        Eigen::Vector3f P_center, Q_center;
        P_center << COM1[0], COM1[1], COM1[2];
        Q_center << COM2[0], COM2[1], COM2[2];
        auto t_found = Q_center - R_found * P_center;

        // start shifting
        for (int i=0;i<pos1.size();i++){
            Eigen::Vector3f vec;
            vec << pos1[i][0], pos1[i][1], pos1[i][2];

            vec = R_found * vec + t_found;

            pos1[i][0] = vec(0);
            pos1[i][1] = vec(1);
            pos1[i][2] = vec(2);
        }

        iteration++;

        if (iteration > 100){
            std::cout << "did not converge." << "\n";
            break;
        }
    }
}

void MeshActions::ChangeMeshWindingOrder(CommandLineArguments& cmd){
    std::string inputfname, outputfname;

    cmd.readString("i", CommandLineArguments::Keys::Required, inputfname);
    cmd.readString("o", CommandLineArguments::Keys::Optional, outputfname);

    Mesh m;
    MeshTools::readPLYlibr(inputfname, m);
    MeshTools::ChangeWindingOrder(m);
    m.CalcVertexNormals();

    MeshTools::writePLY(outputfname, m);
}

void MeshActions::FindFaceNormals(CommandLineArguments& cmd){
    std::string inputfname, outputfname="facenormals.out";

    cmd.readString("i", CommandLineArguments::Keys::Required, inputfname);
    cmd.readString("o", CommandLineArguments::Keys::Optional, outputfname);

    // calculate face normals
    Mesh m;
    MeshTools::readPLYlibr(inputfname, m);
    
    // normals and ares
    std::vector<Real> areas;
    std::vector<Real3> faceNormals;
    MeshTools::CalculateTriangleAreasAndFaceNormals(m, areas, faceNormals);

    std::ofstream ofs;
    ofs.open(outputfname);
    for (int i=0;i<faceNormals.size();i++){
        for (int j=0;j<3;j++){
            ofs << faceNormals[i][j] << " ";
        }
        ofs << "\n";
    }
    ofs.close();
}