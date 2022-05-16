#include "Curvature.h"

Curvature::Curvature(CurvatureInput& input)
:pack_(input.pack)
{
    input.pack.ReadVectorString("outputs", ParameterPack::KeyType::Optional, OutputNames_);
    input.pack.ReadVectorString("outputNames", ParameterPack::KeyType::Optional, OutputFileNames_);
    input.pack.ReadString("name", ParameterPack::KeyType::Optional, name_);

    outputs_.registerOutputFunc("curvature", [this](std::string name) -> void { this -> printCurvature(name);});
    outputs_.registerOutputFunc("principaldir", [this](std::string name) -> void {this -> printPrincipalDir(name);});
    outputs_.registerOutputFunc("FaceCurvature", [this](std::string name)-> void {this -> printFaceCurvature(name);});

    ASSERT(( OutputNames_.size() == OutputFileNames_.size()), "The number of outputs does not agree with number of output files.");
}

void Curvature::initialize(Mesh& mesh)
{
    const auto& vertices = mesh.getvertices();
    int nv = vertices.size();

    CurvaturePerVertex_.clear();
    CurvaturePerVertex_.resize(nv);

    avgCurvaturePerVertex_.clear();
    avgCurvaturePerVertex_.resize(nv ,0.0);

    GaussCurvaturePerVertex_.clear();
    GaussCurvaturePerVertex_.resize(nv,0.0);

    principalDir1_.clear();
    principalDir1_.resize(nv);

    principalDir2_.clear();
    principalDir2_.resize(nv);
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

void Curvature::printFaceCurvature(std::string name)
{
    std::ofstream ofs;
    ofs.open(name);

    ofs << "# TriangleIndex  Curvature\n";
    for (int i=0;i<FaceCurvature_.size();i++)
    {
        ofs << i+1 << " " << FaceCurvature_[i] << "\n";
    }
    ofs.close();
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

void Curvature::CalculateFaceCurvature(Mesh& mesh, const std::vector<Real>& VertexCurvature, std::vector<Real>& FaceCurvature)
{
    const auto& triangle = mesh.gettriangles();
    const auto& vertices = mesh.getvertices();

    FaceCurvature.clear();
    FaceCurvature.resize(triangle.size(),0.0);

    #pragma omp parallel for
    for (int i=0;i<triangle.size();i++)
    {
        auto& t = triangle[i].triangleindices_;

        Real Face_C = 0.0;
        for (int j=0;j<3;j++)
        {
            Face_C += VertexCurvature[t[j]];
        }

        Face_C = Face_C / 3.0;

        FaceCurvature[i] = Face_C;
    }
}