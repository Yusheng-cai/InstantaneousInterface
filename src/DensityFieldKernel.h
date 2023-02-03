#include "GPUArray.cuh"
#include "tools/CommonTypes.h"
#include <thrust/sort.h>
#include <thrust/execution_policy.h>

namespace DensityKernel
{
    using Real = CommonTypes::Real;
    using real3 = Real[3];
    using ind3   = int[3];

    void CalculateInstantaneousField(GPUArray2d<Real>& vector_field_neighbors, GPUArray2d<Real>& atom_positions, GPUArray2d<int>& vector_field_neighbor_index, \
                                     Real inv_sigma2, Real prefactor, GPUArray1d<Real>& box, GPUArray1d<Real>& dx, GPUArray1d<int>& N, GPUArray3d<Real>& insta_field, GPUArray3d<Real>& field, \
                                     GPUArray2d<int>& neighbor_index, int num_atoms, int num_threads=128);

    void CalculateInstantaneousFieldProperty(GPUArray2d<Real>& vector_field_neighbors, GPUArray2d<Real>& atom_positions, GPUArray1d<Real>& Property, \
                                     GPUArray2d<int>& vector_field_neighbor_index, Real inv_sigma2, Real prefactor, GPUArray1d<Real>& box, GPUArray1d<Real>& dx, \
                                     GPUArray1d<int>& N, GPUArray3d<Real>& insta_field, GPUArray3d<Real>& field, \
                                     GPUArray2d<int>& neighbor_index, int num_atoms, int num_threads=128);

    void FillGPUArray(GPUArray3d<Real>& arr, Real num);
} 