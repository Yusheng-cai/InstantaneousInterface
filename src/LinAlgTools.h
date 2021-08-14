#pragma once
#include <tools/CommonTypes.h>

#include <cmath>

namespace LinAlg3x3
{
    using Real = CommonTypes::Real;
    using Matrix = CommonTypes::Matrix;
    using Real3= CommonTypes::Real3;

    Real3 CrossProduct(const Real3& v1, const Real3& v2);
    Real DotProduct(const Real3& v1, const Real3& v2);
    Real norm(const Real3& v1);
    void normalize(Real3& v1);
    Matrix dyad(const Real3& v1, const Real3& v2);

    // LDL^T decomposition of a symmetric positive definite matrix (and some
    // other symmetric matrices, but fragile since we don't do pivoting).
    // Like Cholesky, but no square roots, which is important for small N.
    // Reads diagonal and upper triangle of matrix A.
    // On output, lower triangle of A holds LD, while rdiag holds D^-1.
    // Algorithm from Golub and van Loan.
    bool ldltdc(Matrix& A, Real3 rdiag);

    // Solve Ax=b after ldltdc.  x is allowed to be the same as b. 
    void ldltsl(const Matrix& A, const Real3& rdiag, const Real3& b, Real3& x); 
}