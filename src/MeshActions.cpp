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
    std::string inputfname, outputfname="CurvatureFlow.ply", iterations, lambdadt, decimate="true", numBoundarySmooth="0";
    Real3 box;
    ParameterPack pack;
    bool pbcOutput = false;

    cmd.readString("i", CommandLineArguments::Keys::Required, inputfname);
    cmd.readString("o", CommandLineArguments::Keys::Optional, outputfname);
    cmd.readString("iteration", CommandLineArguments::Keys::Required, iterations);
    cmd.readString("lambdadt", CommandLineArguments::Keys::Required, lambdadt);
    cmd.readBool("pbcOutput", CommandLineArguments::Keys::Optional, pbcOutput);
    cmd.readString("Decimate", CommandLineArguments::Keys::Optional, decimate);
    cmd.readString("NumBoundarySmoothing", CommandLineArguments::Keys::Optional, numBoundarySmooth);
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
    if (pbcMesh)
    {
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
    std::string ext = StringTools::ReadFileExtension(outputfname);
    if (ext == "ply"){MeshTools::writePLY(outputfname, verts, face);}
    else if (ext == "stl"){MeshTools::writeSTL(outputfname, mesh);}
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
                if (isPBC){
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
    for (int i=0;i<lattice.size();i++)
    {
        if (hitRefMeshHeightTotal[i].size() != 0)
        {
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
            if (value[minIndex] < hitRefMeshHeightMin[i]){
                int fcIndex = faceIndex[minIndex];
                pointCurvature[i] = fc[fcIndex];
            }
            else {
                pointCurvature[i] = temp;
            }
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

void MeshActions::MeshPlaneIntersection(CommandLineArguments& cmd)
{
    using Intersector = MeshPlaneIntersect<Real, int>;

    std::vector<Intersector::Vec3D> vertices;
    std::vector<Intersector::Face> faces;
    Real3 p, o;

    std::string inputfname;
    std::string outputfname="intersect.out";

    cmd.readValue("i", CommandLineArguments::Keys::Required, inputfname);
    cmd.readValue("o", CommandLineArguments::Keys::Optional, outputfname);
    cmd.readArray("plane", CommandLineArguments::Keys::Required, p);
    cmd.readArray("points", CommandLineArguments::Keys::Required, o);

    Mesh mesh;
    MeshTools::readPLYlibr(inputfname, mesh);

    const auto& v = mesh.getvertices();
    const auto& f = mesh.gettriangles();

    // copy the data to mesh plane intersection code 
    for (int i=0;i<v.size();i++){
        vertices.push_back(v[i].position_);
    }

    for (int i=0;i<f.size();i++){
        faces.push_back(f[i].triangleindices_);
    }

    // create the mesh for mesh-plane intersection
    Intersector::Mesh m(vertices, faces);
    Intersector::Plane plane;
    plane.normal = p;
    plane.origin = o;

    auto result = m.Intersect(plane);

    std::ofstream ofs;
    ofs.open(outputfname);
    for (int i=0;i<result[0].points.size();i++){
        for (int j=0;j<3;j++){
            ofs << result[0].points[i][j] << " ";
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

void MeshActions::FindBoundaryVertices(CommandLineArguments& cmd)
{
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

void MeshActions::CurvatureEvolution(CommandLineArguments& cmd)
{
    std::string inputfname, outputfname="evolved.ply", neighbors, stepsize, k0, maxstep="1e5", tolerance;
    Real3 box;

    // initialize curve ptr
    refineptr r;
    ParameterPack EvolutionPack, curvePack;

    // read input output
    cmd.readString("i", CommandLineArguments::Keys::Required, inputfname);
    cmd.readString("o", CommandLineArguments::Keys::Optional, outputfname);

    // use curve fit 
    bool isPBC = cmd.readArray("box", CommandLineArguments::Keys::Optional, box);
    cmd.readString("neighbors", CommandLineArguments::Keys::Required, neighbors);
    cmd.readString("stepsize", CommandLineArguments::Keys::Required, stepsize);
    cmd.readString("k0", CommandLineArguments::Keys::Required, k0);
    cmd.readString("maxstep", CommandLineArguments::Keys::Optional, maxstep);
    cmd.readString("tolerance", CommandLineArguments::Keys::Optional, tolerance);
    
    // fill parameter pack
    curvePack.insert("neighbors", neighbors);
    curvePack.insert("type", "curvefit");
    curvePack.insert("name", "temp");
    EvolutionPack.insert("Curvature", curvePack);
    EvolutionPack.insert("name", "refine");
    EvolutionPack.insert("stepsize", stepsize);
    EvolutionPack.insert("k0", k0);
    EvolutionPack.insert("maxstep", maxstep);
    EvolutionPack.insert("tolerance", tolerance);

    // refine pointer 
    MeshRefineStrategyInput input = {{EvolutionPack}};
    r = refineptr(MeshRefineStrategyFactory::Factory::instance().create("CurvatureEvolution", input));

    // read input mesh 
    Mesh m;
    MeshTools::readPLYlibr(inputfname, m);

    // set box length for mesh
    if (isPBC){m.setBoxLength(box);}

    // refine the mesh
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

    for (int i=0;i<nf;i++)
    {
        if (MeshTools::IsIsolatedFace(m, i, mapEdgeToFace)){
            IsolatedTriangles.push_back(t[i].triangleindices_);
        }
    }

    std::ofstream ofs;
    ofs.open(outputfname);
    for (int i=0;i<IsolatedTriangles.size();i++)
    {
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

        for (int i=0;i<nf;i++)
        {
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
    for (int i=0;i<f.size();i++){
        if (MeshTools::IsPeriodicTriangle(mesh, i)){
            PeriodicFaceIndices.push_back(i);
        }
        else{
            NonPeriodicFaceIndices.push_back(i);
            newTriangles.push_back(f[i]);
        }
    }

    // requires periodic meshes --> replicates meshes in the 2 directions
    int t1= translate_index[0], t2=translate_index[1];
    int nv = v.size();
    std::vector<vertex> newVertices;
    for (int i=0;i<array_shape[0];i++){
        for (int j=0;j<array_shape[1];j++){
            for (int k=0;k<nv;k++){
                vertex newV = v[k];
                newV.position_[t1] += i * box[t1];
                newV.position_[t2] += j * box[t2];
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
            int num = i * array_shape[1] + j;
            for (int faceIdx : PeriodicFaceIndices){
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
                Real3 diff1 = v[OtherIndex1].position_ - v[max].position_;
                Real3 diff2 = v[OtherIndex2].position_ - v[max].position_;
                Real3 diff3 = v[OtherIndex1].position_ - v[OtherIndex2].position_;

                // check if either 1 of the edges is both periodic in t1 and t2
                bool bothPeriodic1 = (std::abs(diff1[t1]) >= (0.5 * box[t1])) && (std::abs(diff1[t2]) >= (0.5 * box[t2]));
                bool bothPeriodic2 = (std::abs(diff2[t1]) >= (0.5 * box[t1])) && (std::abs(diff2[t2]) >= (0.5 * box[t2]));
                bool bothPeriodic3 = (std::abs(diff3[t1]) >= (0.5 * box[t1])) && (std::abs(diff3[t2]) >= (0.5 * box[t2]));

                if ((! bothPeriodic1) && (! bothPeriodic2) && (! bothPeriodic3)){
                    if (std::abs(diff1[t2]) >= (0.5 * box[t2]) || std::abs(diff2[t2]) >= (0.5 * box[t2])){
                        newt.triangleindices_[maxIndex] += num * nv;
                        if (std::abs(diff1[t2]) >= (0.5*box[t2])){newt.triangleindices_[OtherIndex[0]] += (num+1)*nv;}
                        else{newt.triangleindices_[OtherIndex[0]] += num * nv;}
                        if (std::abs(diff2[t2]) >= (0.5*box[t2])){newt.triangleindices_[OtherIndex[1]] += (num+1)*nv;}
                        else{newt.triangleindices_[OtherIndex[1]] += num * nv;}
                        newTriangles.push_back(newt);
                    }
                }
            }
        }
    }

    // do all the periodic faces  --> cols
    for (int i=0;i<array_shape[1];i++){
        for (int j=0;j<array_shape[0]-1;j++){
            int num = j * array_shape[1] + i;
            for (int faceIdx : PeriodicFaceIndices){
                // define the new triangle  --> first copy everything from the old triangle
                triangle newt = f[faceIdx];

                // obtain the old triangle indices --> these are real vertex indices 
                INT3 vertsIndex = f[faceIdx].triangleindices_;

                // obtain positions --> in t2 dimension
                Real3 pos = {{v[vertsIndex[0]].position_[t1], v[vertsIndex[1]].position_[t1], v[vertsIndex[2]].position_[t1]}};
                int maxIndex = Algorithm::argmax(pos);
                std::vector<int> OtherIndex;
                for (int k=0;k<3;k++){
                    if (k != maxIndex){OtherIndex.push_back(k);}
                }
                int OtherIndex1 = vertsIndex[OtherIndex[0]], OtherIndex2 = vertsIndex[OtherIndex[1]], max = vertsIndex[maxIndex];

                // find all the edge length in this triangle 
                Real3 diff1 = v[OtherIndex1].position_ - v[max].position_;
                Real3 diff2 = v[OtherIndex2].position_ - v[max].position_;
                Real3 diff3 = v[OtherIndex1].position_ - v[OtherIndex2].position_;

                // check if either 1 of the edges is both periodic in t1 and t2
                bool bothPeriodic1 = (std::abs(diff1[t1]) >= (0.5 * box[t1])) && (std::abs(diff1[t2]) >= (0.5 * box[t2]));
                bool bothPeriodic2 = (std::abs(diff2[t1]) >= (0.5 * box[t1])) && (std::abs(diff2[t2]) >= (0.5 * box[t2]));
                bool bothPeriodic3 = (std::abs(diff3[t1]) >= (0.5 * box[t1])) && (std::abs(diff3[t2]) >= (0.5 * box[t2]));

                if ((! bothPeriodic1) && (! bothPeriodic2) && (! bothPeriodic3)){
                    if (std::abs(diff1[t1]) >= (0.5 * box[t1]) || std::abs(diff2[t1]) >= (0.5 * box[t1])){
                        newt.triangleindices_[maxIndex] += num * nv;
                        if (std::abs(diff1[t1]) >= (0.5*box[t1])){newt.triangleindices_[OtherIndex[0]] += (num+array_shape[1])*nv;}
                        else{newt.triangleindices_[OtherIndex[0]] += num * nv;}
                        if (std::abs(diff2[t1]) >= (0.5*box[t1])){newt.triangleindices_[OtherIndex[1]] += (num+array_shape[1])*nv;}
                        else{newt.triangleindices_[OtherIndex[1]] += num * nv;}
                        newTriangles.push_back(newt);
                    }
                }
            }
        }
    }

    // finally fill the gap  
    for (int i=0;i<array_shape[1];i++){
        for (int j=0;j<array_shape[0]-1;j++){
            int num = j * array_shape[1] + i;
            for (int faceIdx : PeriodicFaceIndices){
                // define the new triangle  --> first copy everything from the old triangle
                triangle newt = f[faceIdx];

                // obtain the old triangle indices --> these are real vertex indices 
                INT3 vertsIndex = f[faceIdx].triangleindices_;

                // obtain positions --> in t2 dimension
                Real3 pos = {{v[vertsIndex[0]].position_[t1], v[vertsIndex[1]].position_[t1], v[vertsIndex[2]].position_[t1]}};
                int maxIndex = Algorithm::argmax(pos);
                std::vector<int> OtherIndex;
                for (int k=0;k<3;k++){
                    if (k != maxIndex){OtherIndex.push_back(k);}
                }
                int OtherIndex1 = vertsIndex[OtherIndex[0]], OtherIndex2 = vertsIndex[OtherIndex[1]], max = vertsIndex[maxIndex];

                // find all the edge length in this triangle 
                Real3 diff1 = v[OtherIndex1].position_ - v[max].position_;
                Real3 diff2 = v[OtherIndex2].position_ - v[max].position_;
                Real3 diff3 = v[OtherIndex1].position_ - v[OtherIndex2].position_;

                // check if either 1 of the edges is both periodic in t1 and t2
                bool bothPeriodic1 = (std::abs(diff1[t1]) >= (0.5 * box[t1])) && (std::abs(diff1[t2]) >= (0.5 * box[t2]));
                bool bothPeriodic2 = (std::abs(diff2[t1]) >= (0.5 * box[t1])) && (std::abs(diff2[t2]) >= (0.5 * box[t2]));
                bool bothPeriodic3 = (std::abs(diff3[t1]) >= (0.5 * box[t1])) && (std::abs(diff3[t2]) >= (0.5 * box[t2]));

                if (bothPeriodic1 || bothPeriodic2 || bothPeriodic3){
                    if (std::abs(diff1[t1]) >= (0.5 * box[t1]) || std::abs(diff2[t1]) >= (0.5 * box[t1])){
                        newt.triangleindices_[maxIndex] += num * nv;
                        if (std::abs(diff1[t1]) >= (0.5*box[t1])){newt.triangleindices_[OtherIndex[0]] += (num+array_shape[1])*nv;}
                        else{newt.triangleindices_[OtherIndex[0]] += num * nv;}
                        if (std::abs(diff2[t1]) >= (0.5*box[t1])){newt.triangleindices_[OtherIndex[1]] += (num+array_shape[1])*nv;}
                        else{newt.triangleindices_[OtherIndex[1]] += num * nv;}
                        newTriangles.push_back(newt);
                    }
                }
            }
        }
    }
    

    v.clear();
    v.insert(v.end(), newVertices.begin(), newVertices.end());
    f.clear();
    f.insert(f.end(), newTriangles.begin(), newTriangles.end());

    MeshTools::writePLY(outputfname, mesh);
}

void MeshActions::TriangleAngleDistribution(CommandLineArguments& cmd){
    using Range = CommonTypes::Real2;
    using Binptr= std::unique_ptr<Bin>;

    std::string inputfname, outputfname="dist.out";
    Real3 box;
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

    for (int i=0;i<angles.size();i++){
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
    std::string inputfname, outputfname="cleaned.ply";

    Real3 box;
    Real edgeLength;
    int iterations=10;
    bool verbose=false;
    cmd.readValue("i", CommandLineArguments::Keys::Required, inputfname);
    cmd.readValue("iterations", CommandLineArguments::Keys::Optional, iterations);
    cmd.readValue("o", CommandLineArguments::Keys::Optional, outputfname);
    cmd.readBool("verbose", CommandLineArguments::Keys::Optional, verbose);
    bool edgeLengthRead = cmd.readValue("edgeLengthCutoff", CommandLineArguments::Keys::Optional, edgeLength);
    bool isPBC = cmd.readArray("box", CommandLineArguments::Keys::Optional, box);

    Mesh m;
    MeshTools::readPLYlibr(inputfname, m);
    if (isPBC){
        m.setBoxLength(box);
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
        ShortEdgeRemoval edgeRemove(m);
        int numCollapsed = edgeRemove.calculate(edgeLength);

        int facesafter = m.gettriangles().size(); 
        if (verbose){
            std::cout << "Collapsed edges = " << numCollapsed << " facesAfter = " << facesafter << "\n";
        }
        int verticesafter = m.getvertices().size();
        if (verticesafter == verticesbefore){
            break;
        }

        if (count > iterations){
            std::cout << "The clean up process did not converge, but exit before it exceed max iteration number " << iterations << "\n";
            break;
        }

        count ++;
    }

    MeshTools::RemoveIsolatedVertices(m);
    MeshTools::RemoveDuplicatedFaces(m);
    MeshTools::RemoveIsolatedFaces(m);

    MeshTools::writePLY(outputfname, m);
}