/* 
 * Clarity is Copyright 2008 Center for Integrated Systems for Microscopy, 
 * Copyright 2008 University of North Carolina at Chapel Hill.
 *
 * Clarity is free software; you can redistribute it and/or modify it under
 * the terms of the GNU Public License as published by the Free Software 
 * Foundation; either version 2 of the License, or (at your option) any 
 * later version.
 *
 * This program is distributed in the hope that it will be useful, but 
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 * or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License 
 * for more details.
 *
 * You should have received a copy of the GNU General Public License along 
 * with this program; if not, write to the Free Software Foundation, Inc., 
 * 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA. You can also find 
 * the GPL on the GNU web site (http://www.gnu.org/copyleft/gpl.html).
 *
 * File name: ComputePrimitivesGPU.cu
 * Author: Cory Quammen <cquammen@cs.unc.edu>
 */


#include <stdio.h>
#include <cuda.h>
#include <cuda_runtime_api.h>
#include "ComputePrimitivesGPU.h"


#define DEFAULT_BLOCKS 64
#define DEFAULT_THREADS_PER_BLOCK 128

#define CLARITY_REDUCE_BLOCKS_ENV            "CLARITY_REDUCE_BLOCKS"
#define CLARITY_REDUCE_THREADS_PER_BLOCK_ENV "CLARITY_REDUCE_THREADS_PER_BLOCK"

int getReduceBlocks() {
  int numBlocks = DEFAULT_BLOCKS;
  char *blocksString = getenv(CLARITY_REDUCE_BLOCKS_ENV);
  if (blocksString) {
    numBlocks = atoi(blocksString);
  }
  
  return numBlocks;
}


int getReduceThreadsPerBlock() {
  int numThreadsPerBlock = DEFAULT_THREADS_PER_BLOCK;
  char *threadsPerBlockString = getenv(CLARITY_REDUCE_THREADS_PER_BLOCK_ENV);
  if (threadsPerBlockString) {
    numThreadsPerBlock = atoi(threadsPerBlockString);
  }
  
  return numThreadsPerBlock;
}


#define CLARITY_MAP_BLOCKS_ENV            "CLARITY_MAP_BLOCKS"
#define CLARITY_MAP_THREADS_PER_BLOCK_ENV "CLARITY_MAP_THREADS_PER_BLOCK"

int getMapBlocks() {
  int numBlocks = DEFAULT_BLOCKS;
  char *blockString = getenv(CLARITY_MAP_BLOCKS_ENV);
  if (blockString) {
    numBlocks = atoi(blockString);
  }

  return numBlocks;
}


int getMapThreadsPerBlock() {
  int numThreadsPerBlock = DEFAULT_THREADS_PER_BLOCK;
  char *threadsPerBlockString = getenv(CLARITY_MAP_THREADS_PER_BLOCK_ENV);
  if (threadsPerBlockString) {
    numThreadsPerBlock = atoi(threadsPerBlockString);
  }

  return numThreadsPerBlock;
}


__global__
void
ReduceSumKernelGPU(float* blockResults, float* data, int n) {
  
  extern __shared__ float accumulator[];
  int tid  = blockDim.x*blockIdx.x + threadIdx.x;
  int incr = gridDim.x*blockDim.x;

  accumulator[threadIdx.x] = 0.0f;
    
  for (int i = tid; i < n; i += incr) {
    // All reads should be coalesced with this pattern.
    accumulator[threadIdx.x] += data[i];
  }
    
  // Reduce the values in shared memory.
  for (int d = blockDim.x >> 1; d > 0; d >>= 1) {
    __syncthreads(); // Make sure all data is read before moving on.
      
    // No bank conflicts in shared memory here.
    if (threadIdx.x < d)
      accumulator[threadIdx.x] += accumulator[threadIdx.x+d];
  }
  __syncthreads();
    
  // Only thread 0 writes the sum to memory.
  if (threadIdx.x == 0)
    blockResults[blockIdx.x] = accumulator[0];
}


extern "C"
void
Clarity_ReduceSumGPU(float* result, float* buffer, int n) {
  
  // Set up device call configuration.
  dim3 gridSize(getReduceBlocks());
  dim3 blockSize(getReduceThreadsPerBlock());
  size_t sharedSize = sizeof(float)*blockSize.x;
  
  // Allocate memory on the device for block-wise partial 
  // reductions computed by the kernel.
  float *blockResultsDev = NULL;
  cudaMalloc((void**)&blockResultsDev, sizeof(float)*gridSize.x);
  
  ReduceSumKernelGPU<<<gridSize, blockSize, sharedSize>>>
    (blockResultsDev, buffer, n);
  
  // Read the partial sums from the blocks back to the host.
  float* blockResultsHost = (float*) malloc(sizeof(float)*gridSize.x);
  cudaMemcpy(blockResultsHost, blockResultsDev, 
	     sizeof(float)*gridSize.x, cudaMemcpyDeviceToHost);
  
  // Add up the results
  *result = 0.0f;
  for (int i = 0; i < gridSize.x; i++) {
    *result += blockResultsHost[i];
  }
  
  free(blockResultsHost);
  cudaFree(blockResultsDev);
}


__global__
void
MultiplyArraysComponentWiseKernelGPU(float* result, float* a, float* b, int n) {
  int tid  = blockDim.x*blockIdx.x + threadIdx.x;
  int incr = gridDim.x*blockDim.x;
  
  for (int i = tid; i < n; i += incr) {
    result[i] = a[i] * b[i];
  }
}


void
Clarity_MultiplyArraysComponentWiseGPU(float* result, float* a, float* b, int n) {
  
  // Set up device call configuration.
  dim3 gridSize(getMapBlocks());
  dim3 blockSize(getMapThreadsPerBlock());

  MultiplyArraysComponentWiseKernelGPU<<<gridSize, blockSize>>>
    (result, a, b, n);
}


__global__
void
DivideArraysComponentWiseKernelGPU(float* result, float* a, float* b, float value, int n) {
  
  int tid  = blockDim.x*blockIdx.x + threadIdx.x;
  int incr = gridDim.x*blockDim.x;
  
  for (int i = tid; i < n; i += incr) {
    if (fabs(b[i]) < 1e-5) {
      result[i] = value;
    } else {
      result[i] = a[i] / b[i];
    }
  }
}


void
Clarity_DivideArraysComponentWiseGPU(float* result, float* a, float* b, float value, int n) {
  
  // Set up device call configuration.
  dim3 gridSize(getMapBlocks());
  dim3 blockSize(getMapThreadsPerBlock());
  
  DivideArraysComponentWiseKernelGPU<<<gridSize, blockSize>>>
    (result, a, b, value, n);

  cudaError error = cudaThreadSynchronize();
  if (error != cudaSuccess) {
    fprintf(stderr, "CUDA error: %s in file '%s' in line %i : %s.\n",
            "Clarity_DivideArraysComponentWiseGPU failed", __FILE__, __LINE__,
            cudaGetErrorString(error));
  }

}


__global__
void
ScaleArrayKernelGPU(float* result, float* a, int n, float scale) {

  int tid  = blockDim.x*blockIdx.x + threadIdx.x;
  int incr = gridDim.x*blockDim.x;
  
  for (int i = tid; i < n; i += incr) {
    result[i] = a[i] * scale;
  }
}


extern "C"
void
Clarity_ScaleArrayGPU(float* result, float* a, int n, float scale) {
  
  // Set up device call configuration.
  dim3 gridSize(getMapBlocks());
  dim3 blockSize(getMapThreadsPerBlock());
  
  ScaleArrayKernelGPU<<<gridSize, blockSize>>>
    (result, a, n, scale);

  cudaError error = cudaThreadSynchronize();
  if (error != cudaSuccess) {
    fprintf(stderr, "CUDA error: %s in file '%s' in line %i : %s.\n",
            "Clarity_ScaleArrayGPU failed", __FILE__, __LINE__,
            cudaGetErrorString(error));
  }
}
