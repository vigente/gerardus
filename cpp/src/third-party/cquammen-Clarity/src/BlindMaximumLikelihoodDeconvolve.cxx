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
 * File name: BlindMaximumLikelihoodDeconvolve.cxx
 * Author: Cory Quammen <cquammen@cs.unc.edu>
 */


#include "Clarity.h"

#include <cstdio>
#include <cstdlib>

#ifdef BUILD_WITH_OPENMP
#include <omp.h>
#endif // BUILD_WITH_OPENMP

#include "ComputePrimitives.h"
#include "Convolve.h"
#include "FFT.h"
#include "MaximumLikelihoodDeconvolve.h"
#include "Memory.h"

extern bool g_CUDACapable;

#ifdef TIME
#include <iostream>
#include "Stopwatch.h"

static Stopwatch totalTimer("BlindMaximumLikelihood filter (total time)");
static Stopwatch transferTimer("BlindMaximumLikelihood filter (transfer time)");
#endif

ClarityResult_t
Clarity_BlindMaximumLikelihoodUpdate(
  int nx, int ny, int nz, float* inImage, float energy, 
  float* o_k, float* h_k, float* o_k_FT, float* h_k_FT,
  float* tmp, float* tmp1_FT, float* tmp2_FT, float* outImage) {

  ClarityResult_t result = CLARITY_SUCCESS;

  // Normalize h_k to preserve energies.
  result = Clarity_NormalizeArray(h_k, h_k, nx*ny*nz);

  // Compute Fourier transform of image guess.
  result = Clarity_FFT_R2C_float(nx, ny, nz, o_k, o_k_FT);

  // Compute Fourier transform of kernel guess.
  result = Clarity_FFT_R2C_float(nx, ny, nz, h_k, h_k_FT);

  // Modulate Fourier transforms of image guess and kernel guess.
  Clarity_Modulate(nx, ny, nz, o_k_FT, h_k_FT, tmp1_FT);

  // Take inverse Fourier transform
  result = Clarity_FFT_C2R_float(nx, ny, nz, tmp1_FT, tmp);

  // Divide input image by convolution of current guess with kernel guess.
  result = Clarity_DivideArraysComponentWise(tmp, inImage, tmp, 0.0f,
					     nx*ny*nz);

  // Take Fourier transform of ratio image.
  result = Clarity_FFT_R2C_float(nx, ny, nz, tmp, tmp1_FT);

  // Convolve the ratio image with h_k.
  Clarity_Modulate(nx, ny, nz, tmp1_FT, h_k_FT, tmp2_FT);

  // Take inverse transform of result.
  result = Clarity_FFT_C2R_float(nx, ny, nz, tmp2_FT, tmp);

  // Multiply o_k with convolution result.
  result = Clarity_MultiplyArraysComponentWise(outImage, o_k, tmp, nx*ny*nz);

  // Convolve the ratio image with o_k.
  Clarity_Modulate(nx, ny, nz, tmp1_FT, o_k_FT, tmp2_FT);

  // Take inverse transform of result.
  result = Clarity_FFT_C2R_float(nx, ny, nz, tmp2_FT, tmp);

  // Multiply h_k with convolution result.
  //result = Clarity_MultiplyArraysComponentWise(h_k, h_k, tmp, nx*ny*nz);

  return result;
}


ClarityResult_t 
Clarity_BlindMaximumLikelihoodDeconvolveCPU(
  float* outImage, float* inImage, float* kernelImage,
  int nx, int ny, int nz, unsigned iterations) {

  ClarityResult_t result = CLARITY_SUCCESS;

  // Yikes, lots of RAM used here.
  float *o_k, *h_k, *o_k_FT, *h_k_FT, *tmp, *tmp1_FT, *tmp2_FT;

  // Set up the array holding the current image guess.
  result = Clarity_Real_Malloc((void**)&o_k, sizeof(float), nx, ny, nz);

  // We'll just alias the kernel guess. This means the kernel image will
  // be clobbered, but that is okay because it is in a padded buffer already.
  h_k = kernelImage;

  // Set up array holding Fourier transform of o_k.
  result = Clarity_Complex_Malloc((void**)&o_k_FT, sizeof(float), nx, ny, nz);

  // Set up array holding Fourier transform of h_k.
  result = Clarity_Complex_Malloc((void**)&h_k_FT, sizeof(float), nx, ny, nz);

  // Storage for intermediate arrays.
  result = Clarity_Real_Malloc((void**) &tmp, sizeof(float), nx, ny, nz);
  result = Clarity_Complex_Malloc((void**)&tmp1_FT, sizeof(float), nx, ny, nz);
  result = Clarity_Complex_Malloc((void**)&tmp2_FT, sizeof(float), nx, ny, nz);

  // Compute original energy in the image
  float energy;
  Clarity_ReduceSum(&energy, inImage, nx*ny*nz);

  // Iterate
  for (unsigned k = 0; k < iterations; k++) {
    float* o_k_cur = (k == 0 ? inImage : o_k);
    float* o_k_new = (k == iterations-1 ? outImage : o_k);

    // Update the image.
    result = Clarity_BlindMaximumLikelihoodUpdate(nx, ny, nz, inImage, energy,
      o_k_cur, h_k, o_k_FT, h_k_FT, tmp, tmp1_FT, tmp2_FT, o_k_new);
  }

  Clarity_Free(o_k); Clarity_Free(o_k_FT);  Clarity_Free(h_k_FT);
  Clarity_Free(tmp); Clarity_Free(tmp1_FT); Clarity_Free(tmp2_FT);

  return result;
}

#ifdef BUILD_WITH_CUDA

#include "MaximumLikelihoodDeconvolveGPU.h"

ClarityResult_t 
Clarity_BlindMaximumLikelihoodDeconvolveGPU(
   float* outImage, float* inImage, float* psfImage, 
   int nx, int ny, int nz, unsigned iterations) {

   ClarityResult_t result = CLARITY_SUCCESS;

   // Copy over PSF and take its Fourier transform.
   float* psf = NULL;
   result = Clarity_Real_MallocCopy((void**) &psf, sizeof(float), 
      nx, ny, nz, psfImage);
   if (result != CLARITY_SUCCESS) {
      return result;
   }
   float* psfFT = NULL;
   result = Clarity_Complex_Malloc((void**) &psfFT, sizeof(float), 
      nx, ny, nz);
   if (result != CLARITY_SUCCESS) {
      Clarity_Free(psf);
      return result;
   }
   result = Clarity_FFT_R2C_float(nx, ny, nz, psf, psfFT);
   if (result != CLARITY_SUCCESS) {
      Clarity_Free(psf); Clarity_Free(psfFT);
      return result;
   }
   Clarity_Free(psf);

   // Copy over image.
   float* in = NULL;
   result = Clarity_Real_MallocCopy((void**) &in, sizeof(float), 
      nx, ny, nz, inImage);
   if (result != CLARITY_SUCCESS) {
      Clarity_Free(psfFT);
      return result;
   }

   // Set up the array holding the current guess.
   float* iPtr = NULL;
   result = Clarity_Real_Malloc((void**)&iPtr, sizeof(float), 
      nx, ny, nz);
   if (result != CLARITY_SUCCESS) {
      Clarity_Free(psfFT); Clarity_Free(in);
      return result;
   }

   // Storage for intermediate arrays
   float* s1 = NULL;
   result = Clarity_Real_Malloc((void**)&s1, sizeof(float), 
      nx, ny, nz);
   if (result != CLARITY_SUCCESS) {
      Clarity_Free(psfFT); Clarity_Free(in); Clarity_Free(iPtr);
      return result;
   }
   float* s2 = NULL;
   result = Clarity_Real_Malloc((void**)&s2, sizeof(float), 
      nx, ny, nz);
   if (result != CLARITY_SUCCESS) {
      Clarity_Free(psfFT); Clarity_Free(in); Clarity_Free(iPtr); 
      Clarity_Free(s1);
      return result;
   }

   // Compute original energy in the image
   float energy;
   Clarity_ReduceSum(&energy, in, nx*ny*nz);

   // Iterate
   for (unsigned k = 0; k < iterations; k++) {
      float* currentGuess = (k == 0 ? in : iPtr);
      float* newGuess     = iPtr;

      result = Clarity_MaximumLikelihoodUpdate(nx, ny, nz, 
         in, energy, currentGuess, psfFT, s1, s2, newGuess);
      if (result != CLARITY_SUCCESS) {
         Clarity_Free(psfFT); Clarity_Free(in); 
         Clarity_Free(iPtr); 
         Clarity_Free(s1); Clarity_Free(s2);
         return result;
      }
   }

   // Copy result from device.
   result = Clarity_CopyFromDevice(nx, ny, nz, sizeof(float), 
      outImage, iPtr);

   Clarity_Free(psfFT); Clarity_Free(in); Clarity_Free(iPtr);
   Clarity_Free(s1);    Clarity_Free(s2);

   return result;
}

#endif // BUILD_WITH_CUDA


ClarityResult_t 
Clarity_BlindMaximumLikelihoodDeconvolve(
  float* inImage, Clarity_Dim3 imageDim, 
  float* kernelImage, Clarity_Dim3 kernelDim,
  float* outImage, int iterations) {

  ClarityResult_t result = CLARITY_SUCCESS;

#ifdef TIME
  totalTimer.Start();
#endif

  // Compute working dimensions. The working dimensions are the sum of the
  // image and kernel dimensions. This handles the cyclic nature of convolution
  // using multiplication in the Fourier domain.
  Clarity_Dim3 workDim;
  workDim.x = imageDim.x + kernelDim.x;
  workDim.y = imageDim.y + kernelDim.y;
  workDim.z = imageDim.z + kernelDim.z;
  int workVoxels = workDim.x*workDim.y*workDim.z;

  // Pad the input image to the working dimensions
  float *inImagePad = (float *) malloc(sizeof(float)*workVoxels);
  int zeroShift[] = {0, 0, 0};
  float fillValue = 0.0f;
  Clarity_ImagePadSpatialShift(inImagePad, workDim, inImage, imageDim,
			       zeroShift, fillValue);

  // Pad the kernel to the working dimensions and shift it so that the
  // center of the kernel is at the origin.
  float *kernelImagePad = (float *) malloc(sizeof(float)*workVoxels);
  int kernelShift[] = {-kernelDim.x/2, -kernelDim.y/2, -kernelDim.z/2};
  Clarity_ImagePadSpatialShift(kernelImagePad, workDim, kernelImage, kernelDim,
			       kernelShift, fillValue);

  // Allocate output array
  float *outImagePad = (float *) malloc(sizeof(float)*workVoxels);

#ifdef BUILD_WITH_CUDA
  if (g_CUDACapable) {
    result = Clarity_BlindMaximumLikelihoodDeconvolveGPU(
        outImagePad, inImagePad, kernelImagePad,
	workDim.x, workDim.y, workDim.z, iterations);
  } else
#endif // BUILD_WITH_CUDA
  {
    result = Clarity_BlindMaximumLikelihoodDeconvolveCPU(
        outImagePad, inImagePad, kernelImagePad,
	workDim.x, workDim.y, workDim.z, iterations);
  }

  // Clip the image to the original dimensions.
  Clarity_ImageClip(outImage, imageDim, outImagePad, workDim);
   
  // Free up memory.
  free(inImagePad);
  free(kernelImagePad);
  free(outImagePad);

#ifdef TIME
  totalTimer.Stop();
  std::cout << totalTimer << std::endl;
  std::cout << transferTimer << std::endl;
  totalTimer.Reset();
  transferTimer.Reset();
#endif

  return result;
}
