#include "GPUArray.cuh"
#include <cuda_runtime.h>
#include "DensityFieldKernel.h"

namespace DensityKernel
{
    __device__ void PBCDistance(real3& vec1, real3& vec2, real3& box, real3& dist, Real& dist_sq)
    {
        dist_sq=0.0;
        for (int i=0;i<3;i++){
            Real diff = vec1[i] - vec2[i];
            if (diff < (-0.5 * box[i])) {diff = diff + box[i];}
            else if (diff > (0.5  * box[i]))  {diff -= box[i];}
            dist_sq += (diff * diff);
            dist[i] = diff;
        }
    }

    __device__ Real GaussianKernel(Real rsq, Real inv_sigma2, Real prefactor)
    {
        return prefactor * expf(rsq * inv_sigma2);
    }

    __device__ void AddIndex(ind3& ind1, ind3& ind2, ind3& maxN, ind3& ret){
        for (int i=0;i<3;i++){
            ret[i] = (ind1[i] + ind2[i]) % maxN[i];

            if (ret[i] < 0){
                ret[i] += maxN[i];
            }
        }
    }

    // function on the device that converts positional index to lattice position
    __device__ void ConvertIndexToLatticePos(ind3& ind, real3& dx, real3& pos){
        for (int i=0;i<3;i++){
            pos[i] = ind[i] * dx[i];
        }
    }

    // function on the device that converts a 3d index into a 1d index
    __device__ void Convert3dIndexTo1d(ind3& index_3d, ind3& max_N, int& index_1d)
    {
        index_1d = index_3d[2] * max_N[0] * max_N[1] + index_3d[1] * max_N[0] + index_3d[0];
    }

    // function that sums the 2 vectors and stores it in C 
    __global__ void Sum(Real* A, Real* B, Real* C, int total_size){
        int idx = blockDim.x * blockIdx.x + threadIdx.x;

        if (idx < total_size){
            C[idx] = A[idx] + B[idx];
        }
    }

    __global__ void SelectiveSum(Real* A, Real* B, Real* C, int* indices, int total_size){
        int idx = blockDim.x * blockIdx.x + threadIdx.x;

        if (idx < total_size){
            int index = indices[idx];

            C[index] = A[index] + B[idx];
        }
    }

    __global__ void CalculateInstantaneousFieldKernel(Real* vector_field, Real* atom_positions, int* vector_field_index, \
                                                Real inv_sigma2, Real prefactor, Real* dx, Real* box, int* total_size, int* neighbors_index, int total_neighbors, int total_data)
    {
        // find the index of the current thread 
        int idx = blockDim.x * blockIdx.x + threadIdx.x;

        // if the index is smaller than the current data 
        if (idx < total_data){
            int atom_idx = idx / total_neighbors;
            int neighbor_idx = idx % total_neighbors;

            // copy over all the data
            ind3 pos_index;
            real3 pos_copy;
            ind3 total_size_copy;
            ind3 offset;
            real3 box_copy;
            real3 dx_copy;
            ind3 actual_index;
            real3 lattice_pos;
            real3 r_vec;
            Real dist_sq=0.0;
            int v_index;

            for (int i=0;i<3;i++){
                pos_copy[i] = atom_positions[atom_idx * 3 + i];
                pos_index[i] = floor(pos_copy[i] / dx[i]);
                box_copy[i] = box[i];
                offset[i] = neighbors_index[neighbor_idx * 3 + i];
                total_size_copy[i] = total_size[i];
                dx_copy[i] = dx[i];
            }

            // add index to 
            AddIndex(pos_index, offset, total_size_copy, actual_index);

            // convert 3d index to 1d 
            Convert3dIndexTo1d(actual_index, total_size_copy, v_index);

            // store the vector index 
            vector_field_index[idx] = v_index;

            // convert index to lattice position
            ConvertIndexToLatticePos(actual_index, dx_copy, lattice_pos);

            // calculate the periodic boundary condition distance 
            PBCDistance(pos_copy, lattice_pos, box_copy, r_vec, dist_sq);

            // calculate the gaussian weight and store into vector_field 
            float gaussWeight = GaussianKernel(dist_sq, inv_sigma2, prefactor);
            vector_field[idx] =  gaussWeight;
        }
        else{
            vector_field[idx] = 0.0;
            vector_field_index[idx] = 0;
        }
    }

}

void DensityKernel::CalculateInstantaneousField(GPUArray2d<Real>& vector_field_neighbors, GPUArray2d<Real>& atom_positions, GPUArray2d<int>& vector_field_neighbor_index, \
                                Real inv_sigma2, Real prefactor, GPUArray1d<Real>& box, GPUArray1d<Real>& dx, GPUArray1d<int>& N, \
                                    GPUArray3d<Real>& insta_field, GPUArray3d<Real>& field, GPUArray2d<int>& neighbor_index, int num_atoms, int num_thread)
{
    int NeighborSize = neighbor_index.getSize()[0];
    int TotalSize = NeighborSize * num_atoms;
    int numBlocks=(TotalSize + num_thread) / num_thread;

    // calculate the instantaneous field 
    thrust::fill(vector_field_neighbors.device_vector().begin(), vector_field_neighbors.device_vector().end(), 0.0);

    DensityKernel::CalculateInstantaneousFieldKernel<<<numBlocks, num_thread>>>(vector_field_neighbors.device_data(), atom_positions.device_data(), vector_field_neighbor_index.device_data(), \
                                                                inv_sigma2, prefactor, dx.device_data(), box.device_data(), N.device_data(), neighbor_index.device_data(), \
                                                                NeighborSize, TotalSize);
    thrust::sort_by_key(thrust::device, vector_field_neighbor_index.device_vector().begin(), vector_field_neighbor_index.device_vector().end(), \
                                        vector_field_neighbors.device_vector().begin());

    auto new_end = thrust::reduce_by_key(thrust::device, vector_field_neighbor_index.device_vector().begin(), vector_field_neighbor_index.device_vector().end(), \
                                        vector_field_neighbors.device_vector().begin(), vector_field_neighbor_index.device_vector().begin(), insta_field.device_vector().begin());
    int new_size = new_end.first - vector_field_neighbor_index.device_vector().begin();

    int num_field_Blocks = (new_size + num_thread) / num_thread;
    DensityKernel::SelectiveSum<<<num_field_Blocks, num_thread>>>(field.device_data(), insta_field.device_data(), field.device_data(), vector_field_neighbor_index.device_data(), new_size);
};

void DensityKernel::FillGPUArray(GPUArray3d<Real>& arr, Real num){
    thrust::fill(arr.device_vector().begin(), arr.device_vector().end(), num);
}