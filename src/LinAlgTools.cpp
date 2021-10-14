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

bool LinAlg3x3::ldltdc(Matrix& A, Real3& rdiag)
{
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

// Solve Ax=b after ldltdc.  x is allowed to be the same as b.
void LinAlg3x3::ldltsl(const Matrix& A,const Real3& rdiag,const Real3& b,Real3& x)
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

LinAlg3x3::Matrix LinAlg3x3::GetRotationMatrix(const Real3& v1, const Real3& v2)
{
    Real3 crossProduct = LinAlg3x3::CrossProduct(v1, v2);
    Real norm = LinAlg3x3::norm(crossProduct);
    Real cosine = LinAlg3x3::DotProduct(v1, v2);

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