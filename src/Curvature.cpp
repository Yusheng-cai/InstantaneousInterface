#include "Curvature.h"

Curvature::Curvature(CurvatureInput& input)
:mesh_(input.mesh_)
{
    input.pack.ReadVectorString("outputs", ParameterPack::KeyType::Optional, OutputNames_);
    input.pack.ReadVectorString("outputNames", ParameterPack::KeyType::Optional, OutputFileNames_);
    outputs_.registerOutputFunc("curvature", [this](std::string name) -> void { this -> printCurvature(name);});
    outputs_.registerOutputFunc("ply", [this](std::string name) -> void {this -> printPLYlibr(name);});
    outputs_.registerOutputFunc("principaldir", [this](std::string name) -> void {this -> printPrincipalDir(name);});

    ASSERT(( OutputNames_.size() == OutputFileNames_.size()), "The number of outputs does not agree with number of output files.");

    const auto& vertices = mesh_.getvertices();
    CurvaturePerVertex_.resize(vertices.size());
    avgCurvaturePerVertex_.resize(vertices.size() ,0.0);
    GaussCurvaturePerVertex_.resize(vertices.size(),0.0);

    principalDir1_.resize(vertices.size());
    principalDir2_.resize(vertices.size());
}

void Curvature::printOutput()
{
    for (int i=0;i<OutputNames_.size();i++)
    {
        outputs_.getOutputFuncByName(OutputNames_[i])(OutputFileNames_[i]);
    }
}

void Curvature::printCurvature(std::string name)
{
    std::ofstream ofs_;
    ofs_.open(name);

    ofs_ << "# k1 k2 Average Gaussian" << "\n";
    for(int i=0;i<CurvaturePerVertex_.size();i++)
    {
        ofs_ << CurvaturePerVertex_[i][0] << " " << CurvaturePerVertex_[i][1] << " " << avgCurvaturePerVertex_[i] << " " \
        << GaussCurvaturePerVertex_[i] << "\n";
    }
    ofs_.close();
}

void Curvature::printPLYlibr(std::string name)
{
    const auto& vertices = mesh_.getvertices();
    const auto& triangles = mesh_.gettriangles();

    std::vector<std::array<double,3>> positions(vertices.size());
    std::vector<std::vector<size_t>> fInd(triangles.size());
    std::vector<std::vector<double>> normals(3, std::vector<double>(vertices.size(),0.0));

    std::vector<std::string> directioNames = {"nx", "ny", "nz"};

    for (int i=0;i<vertices.size();i++)
    {
        for (int j=0;j<3;j++)
        {
            positions[i][j] = vertices[i].position_[j];
            normals[j][i] = vertices[i].normals_[j];
        }
    }

    for (int i=0;i<triangles.size();i++)
    {
        for (int j=0;j<3;j++)
        {
            fInd[i].push_back(triangles[i].triangleindices_[j]);
        }
    }

    happly::PLYData plyOut;

    plyOut.addVertexPositions(positions);
    plyOut.addFaceIndices(fInd);

    for (int i=0;i<3;i++)
    {
        plyOut.getElement("vertex").addProperty(directioNames[i], normals[i]);
    }

    // add average curvature information
    plyOut.getElement("vertex").addProperty("acurv", avgCurvaturePerVertex_);
    plyOut.getElement("vertex").addProperty("gcurve", GaussCurvaturePerVertex_);

    plyOut.write(name, happly::DataFormat::ASCII);
}

void Curvature::printPrincipalDir(std::string name)
{
    std::ofstream ofs_;
    ofs_.open(name);

    ofs_ << "#p1x p1y p1z p2x p2y p2z\n";

    for (int i=0;i<principalDir1_.size();i++)
    {
        for (int j=0;j<3;j++)
        {
            ofs_ << principalDir1_[i][j] << " ";
        }
        for (int j=0;j<3;j++)
        {
            if (j<2)
            {
                ofs_ << principalDir2_[i][j] << " ";
            }
            else
            {
                ofs_ << principalDir2_[i][j] << "\n";
            }
        }
    }

    ofs_.close();
}