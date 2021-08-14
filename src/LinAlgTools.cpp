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

bool LinAlg3x3::ldltdc(Matrix& A, Real3 rdiag)
{
    // Special case for small N
    Real d0 = A[0][0];
    rdiag[0] = 1 / d0;

    A[1][0] = A[0][1];
    Real l10 = rdiag[0] * A[1][0];
    Real d1 = A[1][1] - l10 * A[1][0];
    rdiag[1] = 1 / d1;

    Real d2 = A[2][2] - rdiag[0] * std::pow(A[2][0],2.0) - rdiag[1] * std::pow(A[2][1],2.0);
    rdiag[2] = 1 / d2;
    A[2][0] = A[0][2];
    A[2][1] = A[1][2] - l10 * A[2][0];

    return (d0 != 0 && d1 != 0 && d2 != 0);
}

void LinAlg3x3::ldltsl(const Matrix& A, const Real3& rdiag, const Real3& b, Real3& x) 
{
    for (int i = 0; i < 3; i++) {
        Real sum = b[i];
        for (int k = 0; k < i; k++)
            sum -= A[i][k] * x[k];
        x[i] = sum * rdiag[i];
    }
    for (int i = 3 - 1; i >= 0; i--) {
        Real sum = 0;
        for (int k = i + 1; k < 3; k++)
            sum += A[k][i] * x[k];
        x[i] -= sum * rdiag[i];
    }
}