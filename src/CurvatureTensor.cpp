#include "CurvatureTensor.h"

namespace CurvatureRegistry
{
   registry<CurvatureTensor> registerCurvatureTensor("tensor"); 
}

CurvatureTensor::CurvatureTensor(CurvatureInput& input)
:Curvature(input)
{
    bool readff2 = input.pack.ReadString("2ffOutput", ParameterPack::KeyType::Optional, FF2ofsFileName_);
    if (readff2)
    {
        FF2ofs_.open(FF2ofsFileName_);
        ASSERT((FF2ofs_.is_open()), "The file with name " << FF2ofsFileName_ << " is not opened.");
    }

    bool readprincipalDir = input.pack.ReadString("principaldirOutput", ParameterPack::KeyType::Optional, PrincipalDirectionFileName_);

    if (readprincipalDir)
    {
        PrincipalDirectionofs_.open(PrincipalDirectionFileName_);

        ASSERT((PrincipalDirectionofs_.is_open()), "The file with name " << PrincipalDirectionFileName_ << " is not opened.");
    }

    bool readCurvatureDir = input.pack.ReadVectorArrayNumber("curvaturedir", ParameterPack::KeyType::Optional, curvatureDir_);

    if (readCurvatureDir)
    {
        for (int i=0;i<curvatureDir_.size();i++)
        {
            LinAlg3x3::normalize(curvatureDir_[i]);
        }

        bool readcurvatureoutput = input.pack.ReadString("curvatureDirOutput", ParameterPack::KeyType::Optional, curvatureDirOutputName_);

        if (! readcurvatureoutput)
        {
            std::cout << "WARNING: You are calculating curvature direction without outputting it to some file." << std::endl;
        }
        else
        {
            curvatureDirOutputofs_.open(curvatureDirOutputName_);

            ASSERT((curvatureDirOutputofs_.is_open()), "The file with name " << curvatureDirOutputName_ << " is not opened.");
        }

    }
}

void CurvatureTensor::calculate()
{
    mesh_.CalcPerVertexDir();
    mesh_.CalcTriangleAreaAndFacetNormals();

    const auto& triangles = mesh_.gettriangles();
    const auto& vertices  = mesh_.getvertices();
    const auto& pervertexdir1 = mesh_.getPerVertexDir1();
    const auto& pervertexdir2 = mesh_.getPerVertexDir2();
    const auto& triangleArea = mesh_.getTriangleArea();


    curvatureTensorPerTriangle_.resize(triangles.size());
    curvatureTensorPerVertex_.resize(vertices.size());
    curvatureVec_.resize(vertices.size());
    TotalAreaPerVertex_.resize(vertices.size());
    CurvaturePerVertex_tot.resize(vertices.size());
    std::fill(TotalAreaPerVertex_.begin(), TotalAreaPerVertex_.end(),0);

    for (int i=0;i<triangles.size();i++)
    {
        auto& t = triangles[i].triangleindices_; 
        std::array<Real3,3> edges;

        for (int j=0;j<3;j++)
        {
            edges[0][j] = vertices[t[2]].position_[j] - vertices[t[1]].position_[j];
            edges[1][j] = vertices[t[0]].position_[j] - vertices[t[2]].position_[j];
            edges[2][j] = vertices[t[1]].position_[j] - vertices[t[0]].position_[j];
        }

        Real3 U = edges[0];
        LinAlg3x3::normalize(U);
        Real3 N = LinAlg3x3::CrossProduct(edges[0], edges[1]);
        LinAlg3x3::normalize(N);
        Real3 V = LinAlg3x3::CrossProduct(N,U);
        LinAlg3x3::normalize(V); 

        // initialize A and B matrix to be solved
        Eigen::Matrix3d A;
        Eigen::Vector3d b;
        b.fill(0);
        A.fill(0);
 
        for (int j=0;j<3;j++)
        {
            Real ejU = LinAlg3x3::DotProduct(edges[j], U);
            Real ejV = LinAlg3x3::DotProduct(edges[j], V); 

            A(0,0) += ejU*ejU;
            A(0,1) += ejU*ejV;
            A(2,2) += ejV*ejV;

            int id1 = j - 1;
            if ( id1 < 0)
            {
                id1 += 3;
            }

            int id2 = j - 2;
            if ( id2 < 0)
            {
                id2 += 3;
            }
            Real3 diffN;

            for (int k =0;k<3;k++)
            {
                diffN[k] = vertices[t[id1]].normals_[k] - vertices[t[id2]].normals_[k];
            }

            Real diffNU = LinAlg3x3::DotProduct(diffN, U);
            Real diffNV = LinAlg3x3::DotProduct(diffN, V);
            
            b[0] += diffNU*ejU;
            b[1] += diffNU*ejV + diffNV*ejU;
            b[2] += diffNV*ejV;
        }

        A(1,1) = A(0,0) + A(2,2);
        A(1,2) = A(0,1);
        A(2,1) = A(0,1);
        A(1,0) = A(0,1);

        Eigen::EigenSolver<Eigen::Matrix2d> eigensolver;
        Eigen::Vector3d soln = A.bdcSvd(Eigen::ComputeThinU|Eigen::ComputeThinV).solve(b);
        Real3 ans;
        Eigen::Matrix2d triangleMat;
        triangleMat(0,0) = soln[0];
        triangleMat(0,1) = soln[1];
        triangleMat(1,0) = soln[1];
        triangleMat(1,1) = soln[2];
        eigensolver.compute(triangleMat);

        Eigen::Vector2d eig = eigensolver.eigenvalues().real();
 
        for (int j=0;j<3;j++)
        {
            ans[j] = soln[j];
        }

        curvatureTensorPerTriangle_[i] = ans;
        
        for (int j=0;j<3;j++)
        {
            int id_ = t[j];
            auto& u = pervertexdir1[id_];
            auto& v = pervertexdir2[id_];

            Real3 newcurv = projectCurvature(u,v, U, V,ans);
            


            Eigen::Matrix2d tempMat;
            Eigen::EigenSolver<Eigen::Matrix2d> eigensolver;
            tempMat(0,0) = newcurv[0];
            tempMat(0,1) = newcurv[1];
            tempMat(1,0) = newcurv[1];
            tempMat(1,1) = newcurv[2];

            eigensolver.compute(tempMat);
            Eigen::Vector2d eig = eigensolver.eigenvalues().real();

            CurvaturePerVertex_tot[id_].push_back(newcurv);
            
            for (int k = 0;k<3;k++)
            {
                curvatureTensorPerVertex_[id_][k] += newcurv[k] * triangleArea[i];
            }
            TotalAreaPerVertex_[id_] += triangleArea[i];
        }
    } 

    calculatePrincipalCurvatures();
    calculateCurvatureInDir();
    // for (int i=0;i<CurvaturePerVertex_tot.size();i++)
    // {
    //     if (vertices[i].position_[0] < 2.0)
    //     {
    //         std::cout << "Printing curvature vertices for the " << i << "th vertex." << std::endl;
    //         for (int j=0;j<CurvaturePerVertex_tot[i].size();j++)
    //         {
    //             for (int k=0;k<3;k++)
    //             {
    //                 std::cout << CurvaturePerVertex_tot[i][j][k] << " ";
    //             }
    //             std::cout << "\n";
    //         }
    //     }
    // }
}

void CurvatureTensor::calculatePrincipalCurvatures()
{
    const auto& vertices  = mesh_.getvertices();

    // Use z as a reference direection
    Real3 referenceDir = {{0,0,1}};
    PrincipalDirections_.resize(vertices.size());

    for (int i=0;i<vertices.size();i++)
    {
        // find rotation matrix that rotates the normal from [0,0,1] to itself 
        auto& normal = vertices[i].normals_;
        Matrix rotationMat = LinAlg3x3::GetRotationMatrix(referenceDir, normal);

        for (int j=0;j<3;j++)
        {
            curvatureTensorPerVertex_[i][j] = curvatureTensorPerVertex_[i][j]/TotalAreaPerVertex_[i];
        }

        Eigen::Matrix2d mat;
        Eigen::EigenSolver<Eigen::Matrix2d> eigensolver;

        mat(0,0) = curvatureTensorPerVertex_[i][0];
        mat(0,1) = curvatureTensorPerVertex_[i][1];
        mat(1,0) = curvatureTensorPerVertex_[i][1];
        mat(1,1) = curvatureTensorPerVertex_[i][2];

        eigensolver.compute(mat);
        Eigen::Vector2d eigenvalues = eigensolver.eigenvalues().real();

        curvatureVec_[i][0] = eigenvalues[0];
        curvatureVec_[i][1] = eigenvalues[1];

        Eigen::Matrix2d eigenvectors = eigensolver.eigenvectors().real();

        Real3 eigvec1;
        eigvec1.fill(0);
        Real3 eigvec2;
        eigvec2.fill(0);

        eigvec1[0] = eigenvectors(0,0);
        eigvec1[1] = eigenvectors(1,0);
        eigvec2[0] = eigenvectors(0,1);
        eigvec2[1] = eigenvectors(1,1);

        std::array<Real3, 2> eigenvectorPair;

        Real3 eigvec1Rot = LinAlg3x3::MatrixDotVector(rotationMat, eigvec1);
        Real3 eigvec2Rot = LinAlg3x3::MatrixDotVector(rotationMat, eigvec2);

        // eigenvectors in real space
        eigenvectorPair[0] = eigvec1Rot;
        eigenvectorPair[1] = eigvec2Rot;

        PrincipalDirections_[i]=eigenvectorPair;
    }
}

void CurvatureTensor::printOutput()
{
    if (ofs_.is_open())
    {
        ofs_ << "# k1 k2" << std::endl;
        for (int i=0;i<curvatureVec_.size();i++)
        {
            for (int j=0;j<2;j++)
            {
                ofs_ << curvatureVec_[i][j] << " ";
            }
            ofs_ << "\n";
        }
        ofs_.close();
    }

    if (PrincipalDirectionofs_.is_open())
    {
        for (int i=0;i<PrincipalDirections_.size();i++)
        {
            for (int j=0;j<3;j++)
            {
                PrincipalDirectionofs_ << PrincipalDirections_[i][0][j] << " ";
            }
            for (int j=0;j<3;j++)
            {
                PrincipalDirectionofs_ << PrincipalDirections_[i][1][j] << " ";
            }
            PrincipalDirectionofs_ << "\n";
        }
        PrincipalDirectionofs_.close();
    }

    if (FF2ofs_.is_open())
    {
        for (int i=0;i<curvatureTensorPerVertex_.size();i++)
        {
            for (int j=0;j<3;j++)
            {
                FF2ofs_ << curvatureTensorPerVertex_[i][j] << " ";
            }
            FF2ofs_ << "\n";
        }
        FF2ofs_.close();
    }

    if (curvatureDirOutputofs_.is_open())
    {
        const auto& vertices = mesh_.getvertices();

        for (int i=0;i<vertices.size();i++)
        {
            for (int j=0;j<curvaturesInDir_.size();j++)
            {
                curvatureDirOutputofs_ << curvaturesInDir_[j][i] << " ";
            }
            curvatureDirOutputofs_ << "\n";
        }
        curvatureDirOutputofs_.close();
    }
}

void CurvatureTensor::calculateCurvatureInDir()
{
    if (curvatureDir_.size() != 0)
    {
        const auto& vertices = mesh_.getvertices();
        const auto& PerVertexDir1 = mesh_.getPerVertexDir1();
        const auto& PerVertexDir2 = mesh_.getPerVertexDir2();
        curvaturesInDir_.resize(curvatureDir_.size(), std::vector<Real>(vertices.size()));

        for (int i=0;i<curvatureDir_.size();i++)
        {
            auto dir = curvatureDir_[i];
            for (int j=0;j<vertices.size();j++)
            {
                auto& tensor = curvatureTensorPerVertex_[j];
                auto& dir1   = PerVertexDir1[j];
                auto& dir2   = PerVertexDir2[j];

                Real s = LinAlg3x3::DotProduct(dir1, dir);
                Real t = LinAlg3x3::DotProduct(dir2, dir);

                Real c = s*s*tensor[0] + 2*s*t*tensor[1] + tensor[2]*t*t;

                curvaturesInDir_[i][j] = c;
            }
        }
    }
}

CurvatureTensor::Real3 CurvatureTensor::projectCurvature(const Real3& oldu, const Real3& oldv, const Real3& refu, const Real3& refv,const Real3& curvature)
{
    Real3 newu;
    Real3 newv;
 

    Real3 oldN = LinAlg3x3::CrossProduct(oldu, oldv);
    LinAlg3x3::normalize(oldN);

    Real3 refN = LinAlg3x3::CrossProduct(refu, refv);
    LinAlg3x3::normalize(refN);

    Matrix rotationMatrix = LinAlg3x3::GetRotationMatrix(oldN, refN);

    newu = LinAlg3x3::MatrixDotVector(rotationMatrix, oldu);
    LinAlg3x3::normalize(newu);
    newv = LinAlg3x3::MatrixDotVector(rotationMatrix, oldv);
    LinAlg3x3::normalize(newv);

	Real u1 = LinAlg3x3::DotProduct(newu, refu);
    Real v1 = LinAlg3x3::DotProduct(newu, refv);
    Real u2 = LinAlg3x3::DotProduct(newv, refu);
    Real v2 = LinAlg3x3::DotProduct(newv, refv);

    Real3 newcurvature_;
	newcurvature_[0]  = curvature[0] * u1*u1 + curvature[1] * (2.0  * u1*v1) + curvature[2] * v1*v1;
	newcurvature_[1]  = curvature[0] * u1*u2 + curvature[1] * (u1*v2 + u2*v1) + curvature[2] * v1*v2;
	newcurvature_[2]  = curvature[0] * u2*u2 + curvature[1] * (2.0  * u2*v2) + curvature[2] * v2*v2;

    return newcurvature_;
}