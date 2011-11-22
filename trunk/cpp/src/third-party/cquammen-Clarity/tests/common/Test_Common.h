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
 * File name: Test_Common.h
 * Author: Cory Quammen <cquammen@cs.unc.edu>
 */


#ifndef __TEST_COMMON_H_
#define __TEST_COMMON_H_

#include <cmath>
#include <cstdlib>
#include <cstdio>

#define IMG_X 128
#define IMG_Y 128
#define IMG_Z 32

#define KERNEL_X 32
#define KERNEL_Y 32
#define KERNEL_Z 32

#ifndef M_PI
#define M_PI 3.14159265
#endif

/** Determines the size of the image used in the tests.
 *
 * Attempts to read environment variables to determine image size. If one of
 * the environment variables does not exist, the default image size is
 * returned. Environment variables affecting this function are:
 * CLARITY_IMAGE_SIZE_X, CLARITY_IMAGE_SIZE_Y, and CLARITY_IMAGE_SIZE_Z.
 */
void getImageSize(Clarity_Dim3 *imageDims) {
  imageDims->x = IMG_X;
  imageDims->y = IMG_Y;
  imageDims->z = IMG_Z;

  char *xDimString = getenv("CLARITY_IMAGE_SIZE_X");
  char *yDimString = getenv("CLARITY_IMAGE_SIZE_Y");
  char *zDimString = getenv("CLARITY_IMAGE_SIZE_Z");
  if (xDimString && yDimString && zDimString) {
    imageDims->x = atoi(xDimString);
    imageDims->y = atoi(yDimString);
    imageDims->z = atoi(zDimString);
  }

}


/** Determines the size of the convolution kernel used in the tests.
 *
 * Attempts to read environment variables to determine kernel size. If one of
 * the environment variables does not exist, the default kernel size is
 * returned. Environment variables affecting this function are:
 * CLARITY_KERNEL_SIZE_X, CLARITY_KERNEL_SIZE_Y, and CLARITY_KERNEL_SIZE_Z.
 */
void getKernelSize(Clarity_Dim3 *kernelDims) {
  kernelDims->x = KERNEL_X;
  kernelDims->y = KERNEL_Y;
  kernelDims->z = KERNEL_Z;

  char *xDimString = getenv("CLARITY_KERNEL_SIZE_X");
  char *yDimString = getenv("CLARITY_KERNEL_SIZE_Y");
  char *zDimString = getenv("CLARITY_KERNEL_SIZE_Z");
  if (xDimString && yDimString && zDimString) {
    kernelDims->x = atoi(xDimString);
    kernelDims->y = atoi(yDimString);
    kernelDims->z = atoi(zDimString);
  }

}


/**
 * Generates image representing the true signal.
 * 
 * @param dim Return parameter for the dimensions of the true signal.
 * @return Image stored in memory allocated by malloc. The caller is
 * responsible for free'ing this memory.
 */
float *
Test_GenerateTrueImage(Clarity_Dim3 *dim) {
  getImageSize(dim);

  float *image = (float *) malloc(sizeof(float)*dim->x*dim->y*dim->z);

  // Initialize with zeros.
  for (int i = 0; i < dim->x*dim->y*dim->z; i++) {
    image[i] = 0.0f;
  }

  for (int iz = dim->z/4; iz < dim->z - (dim->z/4); iz++) {
    for (int iy = dim->y/4; iy < dim->y - (dim->y/4); iy++) {
      for (int ix = dim->x/4; ix < dim->x - (dim->x/4); ix++) {
        image[iz*dim->x*dim->y + iy*dim->x + ix] = 1.0f;
      }
    }
  }

  return image;
}


/**
 * Generates image representing a Gaussian convolution kernel.
 *
 * @param dim   Return parameter for the dimensions of the kernel.
 * @param sigma Standard deviation of blurring kernel.
 * @return Kernel image stored in memory allocated by malloc. The caller is
 * responsible for free'ing this memory.
 */
float *
Test_GenerateGaussianKernel(Clarity_Dim3 *dim, float sigma) {
  getKernelSize(dim);

  float *kernel = (float *) malloc(sizeof(float)*dim->x*dim->y*dim->z);

  float sum = 0.0f;
  float sigma2 = sigma*sigma;
  for (int iz = 0; iz < dim->z; iz++) {
    float fz = static_cast<float>(iz-(dim->z/2));
    for (int iy = 0; iy < dim->y; iy++) {
      float fy = static_cast<float>(iy-(dim->y/2));
      for (int ix = 0; ix < dim->x; ix++) {
        float fx = static_cast<float>(ix-(dim->x/2));
	float value =
          (1.0f / pow(2.0*M_PI*sigma2, 1.5)) *
          exp(-((fx*fx + fy*fy + fz*fz)/(2*sigma2)));
        kernel[(iz*dim->x*dim->y) + (iy*dim->x) + ix] = value;
	sum += value;
      }
    }
  }

  // Normalize the kernel
  float div = 1.0f / sum;
  for (int i = 0; i < dim->x*dim->y*dim->z; i++) {
    kernel[i] *= div;
  }

  return kernel;

}


/**
 * Reports match between a known deconvolution solution and the
 * deconvolved image.
 *
 * @param inputImage       Known deconvolution solution.
 * @param deconvolvedImage Result of deconvolution algorithm.
 * @param imageDims        Dimensions of inputImage and deconvolvedImage.
 */
void
Test_ReportMatch(float* inputImage, float* deconvolvedImage,
		 Clarity_Dim3 imageDims) {
  float inputSum = 0.0f;
  float deconvolvedSum = 0.0f;
  float sum2 = 0.0f;
  for (int iz = 0; iz < imageDims.z; iz++) {
    for (int iy = 0; iy < imageDims.y; iy++) {
      for (int ix = 0; ix < imageDims.x; ix++) {
	int index = iz*imageDims.x*imageDims.y + iy*imageDims.x + ix;
	float diff = deconvolvedImage[index] - inputImage[index];
	sum2 += diff*diff;
	inputSum += inputImage[index];
	deconvolvedSum += deconvolvedImage[index];
      }
    }
  }

  printf("RMS is: %f\n", 
	 sqrt(sum2/static_cast<float>(imageDims.x*imageDims.y*imageDims.z)));
  printf("Difference in total intensity between images: %f\n",
	 inputSum - deconvolvedSum);

}


#endif // __TEST_COMMON_H_
