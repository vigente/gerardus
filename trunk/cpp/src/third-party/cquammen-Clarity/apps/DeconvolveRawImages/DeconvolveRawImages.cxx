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
 * File name: DeconvolveRawImages.cxx
 * Author: Cory Quammen <cquammen@cs.unc.edu>
 */

/*
 * This is a derivative work of Clarity provided as a third-party
 * library in Gerardus
 *
 *  Minor fixes by Ramon Casero <rcasero@gmail.com>
 */

/* This example shows how to read unsigned short images, run a deconvolution
 * routine on them, and write the result. */

#include <cstdio>
#include <cstdlib>

#include <Clarity.h>

#define BUF_SIZE 1024

int main(int argc, char* argv[]) {
  
  if (argc < 10) {
    printf("Usage: DeconvolveRawImages <image filename> <x> <y> <z> <psf filename> <x> <y> <z> <output file name>\n");
    return 1;
  }

  // Parse command-line arguments
  char *imageFileName = argv[1];
  Clarity_Dim3 imageDim;
  imageDim.x = atoi(argv[2]);
  imageDim.y = atoi(argv[3]);
  imageDim.z = atoi(argv[4]);
  char *psfFileName = argv[5];
  Clarity_Dim3 psfDim;
  psfDim.x = atoi(argv[6]);
  psfDim.y = atoi(argv[7]);
  psfDim.z = atoi(argv[8]);
  char *outputFileName = argv[9];

  // Read in image and PSF
  FILE *inputFp = fopen(imageFileName, "rb");
  FILE *psfFp = fopen(psfFileName, "rb");
  FILE *outputFp = fopen(outputFileName, "wb");
  if (!inputFp || !psfFp || !outputFp) {
    if (!inputFp) 
      printf("ERROR: Could not open image file '%s'\n", imageFileName);
    if (!psfFp) 
      printf("ERROR: Could not open PSF file '%s'\n'", psfFileName);
    if (!outputFp) 
      printf("ERROR: Could not open file '%s' for writing\n", outputFileName);

    if (inputFp)  fclose(inputFp);
    if (psfFp)    fclose(psfFp);
    if (outputFp) fclose(outputFp);
    return -1;
  }
  
  // Everything's okay, let's read.
  int inputImageSize = imageDim.x*imageDim.y*imageDim.z;
  int psfImageSize   = psfDim.x*psfDim.y*psfDim.z;
  unsigned short *inputImage  = new unsigned short[inputImageSize];
  unsigned short *psfImage    = new unsigned short[psfImageSize];

  if ((size_t)inputImageSize 
      != fread(inputImage, sizeof(unsigned short), inputImageSize, inputFp)) {
    printf("Error reading input image file\n");
    return 1;
  }
  fclose(inputFp);
  if ((size_t)psfImageSize 
      != fread(psfImage, sizeof(unsigned short), psfImageSize, psfFp)) {
    printf("Error reading PSF file\n");
    return 1;
  }
  fclose(psfFp);
  
  // Cast to floats
  float *inputImageFloat = new float[inputImageSize];
  for (int i = 0; i < inputImageSize; i++) {
    inputImageFloat[i] = static_cast<float>(inputImage[i]);
  }
  delete[] inputImage;

  float *psfImageFloat = new float[psfImageSize];
  for (int i = 0; i < psfImageSize; i++) {
    psfImageFloat[i] = static_cast<float>(psfImage[i]);
  }
  delete[] psfImage;

  // Allocate output buffer.
  float *outputImageFloat = new float[inputImageSize];

  // Run the deconvolution. We'll use the maximum likelihood method.
  Clarity_MaximumLikelihoodDeconvolve(inputImageFloat, imageDim,
                                      psfImageFloat, psfDim,
                                      outputImageFloat, 5);

  // Cast back to unsigned shorts and write.
  unsigned short *outputImage = new unsigned short[inputImageSize];
  for (int i = 0; i < inputImageSize; i++) {
    outputImage[i] = static_cast<unsigned short>(outputImageFloat[i]);
  }

  fwrite(outputImage, sizeof(unsigned short), inputImageSize, outputFp);
  fclose(outputFp);
  delete[] outputImage;

  // Report what we did.
  printf("Input image file name: %s\n", imageFileName);
  printf("Input image dimensions: %d x %d x %d\n",
         imageDim.x, imageDim.y, imageDim.z);
  printf("Input point-spread function file name: %s\n", psfFileName);
  printf("Input point-spread function dimensions: %d x %d x %d\n",
         psfDim.x, psfDim.y, psfDim.z);
  printf("Output file name: %s\n", outputFileName);

  return 0;
}
