#include "CurvatureTensor.h"

namespace CurvatureRegistry
{
   registry<CurvatureTensor> registerCurvatureTensor("tensor"); 
}

CurvatureTensor::CurvatureTensor(CurvatureInput& input)
:Curvature(input)
{}

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
        Real3 V = LinAlg3x3::CrossProduct(U, N);
        LinAlg3x3::normalize(V); 

        // initialize A and B matrix to be solved
        Eigen::Matrix3d A;
        Eigen::Vector3d b;
        b.fill(0);
        A.fill({});
 
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

        Eigen::Vector3d soln = A.bdcSvd(Eigen::ComputeThinU|Eigen::ComputeThinV).solve(b);
        Real3 ans;

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
            for (int k = 0;k<3;k++)
            {
                curvatureTensorPerVertex_[id_][k] += newcurv[k] * triangleArea[i];
            }
            TotalAreaPerVertex_[id_] += triangleArea[i];
        }
    }

    for (int i=0;i<vertices.size();i++)
    {
        for (int j=0;j<3;j++)
        {
            curvatureTensorPerVertex_[i][j] /= TotalAreaPerVertex_[i];
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
}

CurvatureTensor::Matrix CurvatureTensor::getRotationMatrix(const Real3& vec1, const Real3& vec2)
{
    Real3 crossProduct = LinAlg3x3::CrossProduct(vec1, vec2);
    Real norm = LinAlg3x3::norm(crossProduct);
    Real cosine = LinAlg3x3::DotProduct(vec1, vec2);

    Matrix ret;
    Real denom = 1.0 + cosine;
    Real factor;

    if (std::abs(denom) < 1e-7) { factor = 0.0;}
    else{factor = 1.0/(1.0 + cosine);}


    ret[0][0] = 1 + factor*(-crossProduct[2]*crossProduct[2] - crossProduct[1]*crossProduct[1]);
    ret[0][1] = -crossProduct[2] + factor*crossProduct[0]*crossProduct[1];
    ret[0][2] = crossProduct[1] + factor*crossProduct[0]*crossProduct[2];
    ret[1][0] = crossProduct[2] + factor*crossProduct[0]*crossProduct[1];
    ret[1][1] = 1 + factor*(-crossProduct[2]*crossProduct[2] - crossProduct[0]*crossProduct[0]);
    ret[1][2] = -crossProduct[0] + factor*crossProduct[1]*crossProduct[2];
    ret[2][0] = -crossProduct[1] + factor*crossProduct[2]*crossProduct[0];
    ret[2][1] = crossProduct[0] + factor*crossProduct[2]*crossProduct[1];
    ret[2][2] = 1 + factor*(-crossProduct[1]*crossProduct[1] - crossProduct[0]*crossProduct[0]);

    return ret;
}

CurvatureTensor::Real3 CurvatureTensor::projectCurvature(const Real3& oldu, const Real3& oldv, const Real3& refu, const Real3& refv,const Real3& curvature)
{
    Real3 newu;
    Real3 newv;
 

    Real3 oldN = LinAlg3x3::CrossProduct(oldu, oldv);
    LinAlg3x3::normalize(oldN);

    Real3 refN = LinAlg3x3::CrossProduct(refu, refv);
    LinAlg3x3::normalize(refN);

    Matrix rotationMatrix = getRotationMatrix(oldN, refN);

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
	newcurvature_[1] =  curvature[0] * u1*u2 + curvature[1] * (u1*v2 + u2*v1) + curvature[2] * v1*v2;
	newcurvature_[2]  = curvature[0] * u2*u2 + curvature[1] * (2.0  * u2*v2) + curvature[2] * v2*v2;

    return newcurvature_;
}
