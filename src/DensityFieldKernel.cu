#include "GPUArray.h"
#include <cuda_runtime.h>

namespace DensityKernel
{
    using real3 = float[3];
    using ind3   = int[3];

    __device__ void PBCDistance(real3& vec1, real3& vec2, real3& box, real3& dist, float& dist_sq)
    {
        dist_sq=0.0;
        for (int i=0;i<3;i++){
            float diff = vec1[i] - vec2[i];
            if (diff < (-0.5 * box[i])) {diff = diff + box[i];}
            else if (diff > (0.5  * box[i]))  {diff -= box[i];}
            dist_sq += (diff * diff);
            dist[i] = diff;
        }
    }

    __device__ float GaussianKernel(float rsq, float inv_sigma2, float prefactor)
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

    __device__ void ConvertIndexToLatticePos(ind3& ind, real3& dx, real3& pos){
        for (int i=0;i<3;i++){
            pos[i] = ind[i] * dx[i];
        }
    }

    __device__ void LocalIndexToNeighborOffset(int index, ind3& neighbors, ind3& offset){
        int xdim = 2 * neighbors[0] + 1;
        int ydim = 2 * neighbors[1] + 1;
        int xy = xdim * ydim;

        int zindex = index / xy; 
        int yindex = (index - zindex * xy) / xdim;
        int xindex = index - zindex * xy - yindex * xdim;

        offset[0] = xindex - neighbors[0];
        offset[1] = yindex - neighbors[1];
        offset[2] = zindex - neighbors[2];
    }

    __device__ void Convert3dIndexTo1d(ind3& index_3d, ind3& max_N, int& index_1d)
    {
        index_1d = index_3d[2] * max_N[0] * max_N[1] + index_3d[1] * max_N[0] + index_3d[0];
    }

    __global__ void Sum(float* A, float* B, float* C, int total_size){
        int idx = blockDim.x * blockIdx.x + threadIdx.x;

        if (idx < total_size){
            C[idx] = A[idx] + B[idx];
        }
    }

    __global__ void CalculateInstantaneousField(float* vector_field, float* atom_positions, int* vector_field_index, float inv_sigma2, float prefactor, float* dx, float* box, int* total_size, int* neighbors, int total_data)
    {
        int idx = blockDim.x * blockIdx.x + threadIdx.x;

        if (idx < total_data){
            int total_neighbors = (2 * neighbors[0] + 1) * (2 * neighbors[1] + 1) * (2 * neighbors[2] + 1);
            int atom_idx = idx / total_neighbors;
            int neighbor_idx = idx % total_neighbors;

            // // copy over the positions
            ind3 pos_index;
            real3 pos_copy;
            ind3 total_size_copy;
            ind3 neighbor_size_copy;
            real3 box_copy;
            real3 dx_copy;

            ind3 actual_index;
            real3 lattice_pos;
            real3 dist;
            float dist_sq;

            for (int i=0;i<3;i++){
                pos_copy[i] = atom_positions[atom_idx*3+i];
                pos_index[i] = floor(pos_copy[i] / dx[i]);
                box_copy[i] = box[i];
                neighbor_size_copy[i] = neighbors[i];
                total_size_copy[i] = total_size[i];
                dx_copy[i] = dx[i];
            }

            ind3 offset;
            int v_index;
            LocalIndexToNeighborOffset(neighbor_idx, neighbor_size_copy, offset);
            AddIndex(pos_index, offset, total_size_copy, actual_index);
            Convert3dIndexTo1d(actual_index, total_size_copy, v_index);
            vector_field_index[idx] = v_index;
            ConvertIndexToLatticePos(actual_index, dx_copy, lattice_pos);
            PBCDistance(pos_copy, lattice_pos, box_copy, dist, dist_sq);
            float gauss_weight = GaussianKernel(dist_sq, inv_sigma2, prefactor);
            vector_field[idx] = gauss_weight;
        }
        else{
            vector_field[idx] = 0.0;
            vector_field_index[idx] = 0;
        }
    }
}