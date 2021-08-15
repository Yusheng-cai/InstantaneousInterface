#include "LinAlgTools.h"

LinAlg3x3::Real3 LinAlg3x3::CrossProduct(const Real3& v1, const Real3& v2)
{
    Real3 ret;
    ret.fill(0);

    ret[0] = v1[1]*v2[2] - v1[2]*v2[1];
    ret[1] = v1[2]*v2[0] - v1[0]*v2[2];
    ret[2] = v1[0]*v2[1] - v1[1]*v2[0];

    return ret;
}

LinAlg3x3::Real LinAlg3x3::DotProduct(const Real3& v1, const Real3& v2)
{
    Real ret = 0.0;

    for (int i=0;i<3;i++)
    {
        ret += v1[i]*v2[i];
    }

    return ret;
}

LinAlg3x3::Real LinAlg3x3::norm(const Real3& v1)
{
    Real ret = 0.0;

    for (int i=0;i<3;i++)
    {
        ret += v1[i]*v1[i];
    }

    ret = std::sqrt(ret);

    return ret;
}

LinAlg3x3::Matrix LinAlg3x3::dyad(const Real3& v1, const Real3& v2)
{
    Matrix ret;

    for (int i=0;i<3;i++)
    {
        for (int j=0;j<3;j++)
        {
            ret[i][j] = v1[i]*v2[j];
        }
    }

    return ret;
}

void LinAlg3x3::normalize(Real3& v1)
{
    Real norm_ = norm(v1);

    for (int i=0;i<3;i++)
    {
        v1[i] = v1[i]/norm_;
    }
}

LinAlg3x3::Real3 LinAlg3x3::MatrixDotVector(const Matrix& A, const Real3& v1)
{
    Real3 ans;
    ans.fill(0);

    for (int i=0;i<3;i++)
    {
        for(int j=0;j<3;j++)
        {
            ans[i] += A[i][j]*v1[j];
        }
    }

    return ans;
}