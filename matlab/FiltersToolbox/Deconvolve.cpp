/**
 * DECONVOLVE  MEX interface to the Clarity deconvolution library
 *
 * This function provides a Matlab MEX interface to the Clarity
 * Deconvolution Library [1] developed by CISMM, UNC [2].
 *
 * Clarity uses by default one thread per processor. In a 4-processor
 * machine, we have observed typical run times 1/2 to 1/3 of Matlab's native
 * deconvolution algorithms (e.g. deconvwnr, deconvlucy).
 *
 * Clarity provides a CUDA implementation too. To enable it, read the
 * documentation in the wiki:
 *
 * http://code.google.com/p/gerardus/wiki/EnablingCUDA
 *
 * A = DECONVOLVE(ALGO, B, PSF)
 *
 *   The working hypothesis is that we want to recover A from
 *
 *     B = convolution(A, PSF) + noise
 *
 *   ALGO is a string that selects the deconvolution method:
 *
 *     'w': Clarity_WienerDeconvolve. Wiener filter deconvolution. Matlab
 *     equivalent deconvwnr().
 *
 *     'jvc': Clarity_JansenVanCittertDeconvolve. Classic Jansen-van Cittert
 *     formulation for constrained iterative deconvolution.
 *
 *     'ml': Clarity_MaximumLikelihoodDeconvolve. Maximum-likelihood
 *     deconvolution method cited in [3].
 *
 *   B is the observed image, affected by blurring and noise. B is typically
 *   a 2- or 3-dimensional array.
 *
 *   PSF is an array with an estimate of the point spread function (PSF).
 *
 *   A is the deconvolved image. A is an array with the same dimensions as
 *   B.
 *
 * A = DECONVOLVE('w', B, PSF, EPSILON)
 *
 *   EPSILON is a scalar. From the Clarity documentation: "Constant standing
 *   in place of the ratio between power spectra of noise and the power
 *   spectra of the underlying image, which are unknown parameters. In
 *   practice, acts as a smoothing factor. Typically set in the range 0.001
 *   to 0.1". By default, EPSILON=0.01.
 *
 * A = DECONVOLVE('jvc', B, PSF, ITER)
 * A = DECONVOLVE('ml', B, PSF, ITER)
 *
 *   ITER is a scalar. Number of iterations for the algorithm. By default,
 *   ITER=10.
 *
 * [1] http://cismm.cs.unc.edu/resources/software-manuals/clarity-deconvolution-library/
 * [2] Computer Integrated Systems for Microscopy and Manipulation,
 *     University of North Carolina, http://cismm.cs.unc.edu/
 * [3] J.B. Sibarita, Deconvolution microscopy, Adv. Biochem.
 * Engin./Biotechnology (2005) 95: 201-243.
*/

/**
 * Author: Ramon Casero <rcasero@gmail.com>
 * Copyright © 2011 University of Oxford
 * Version: 0.1.0
 * $Rev$
 * $Date$
 * 
 * University of Oxford means the Chancellor, Masters and Scholars of
 * the University of Oxford, having an administrative office at
 * Wellington Square, Oxford OX1 2JD, UK. 
 *
 * This file is part of Gerardus.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details. The offer of this
 * program under the terms of the License is subject to the License
 * being interpreted in accordance with English Law and subject to any
 * action against the University of Oxford being under the jurisdiction
 * of the English Courts.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/* mex headers */
#include <math.h>
#include <matrix.h>
#include <mex.h>

/* C++ headers */
#include <iostream>
#include <cstring>

/* Clarity headers */
#include <Clarity.h>

void mexFunction(int nlhs, mxArray *plhs[], 
		 int nrhs, const mxArray *prhs[]) {

  // check number of input and output arguments
  if (nrhs < 3) {
    mexErrMsgTxt("Too few input arguments");
  }
  if (nrhs > 4) {
    mexErrMsgTxt("Too many input arguments");
  }
  if (nlhs > 1) {
    mexErrMsgTxt("Too many output arguments");
  }

  // get string with deconvolution algorithm
  if (mxGetClassID(prhs[0]) != mxCHAR_CLASS) {
    mexErrMsgTxt("First argument must be a string to select one deconvolution algorithm");
  }
  
  // get type of filter
  char *algorithmName = mxArrayToString(prhs[0]);
  if (algorithmName == NULL) {
    mexErrMsgTxt("Invalid ALGORITHM string");
  }
  
  // get input image class (double, uint8, etc)
  if ((mxGetClassID(prhs[1]) != mxSINGLE_CLASS) ||
      (mxGetClassID(prhs[2]) != mxSINGLE_CLASS)) {
    mexErrMsgTxt("deconvolution library Clarity requires images given as float type (single type in Matlab)");
  }
  
  // get image and PSF dimensions
  Clarity_Dim3 imageDims; // Image dimensions
  Clarity_Dim3 kernelDims; // PSF dimensions
  mwSize imageNDims = mxGetNumberOfDimensions(prhs[1]);
  mwSize kernelNDims = mxGetNumberOfDimensions(prhs[2]);
  const mwSize *imageMatDims = mxGetDimensions(prhs[1]);
  const mwSize *kernelMatDims = mxGetDimensions(prhs[2]);

  imageDims.y = imageMatDims[0]; // row = y-coordinate
  imageDims.x = imageMatDims[1]; // col = x-coordinate
  if (imageNDims < 3) {
    imageDims.z = 1; // 1 slice
  } else {
    imageDims.z = imageMatDims[2]; // slice = z-coordinate
  }

  kernelDims.y = kernelMatDims[0]; // row = y-coordinate
  kernelDims.x = kernelMatDims[1]; // col = x-coordinate
  if (kernelNDims < 3) {
    kernelDims.z = 1; // 1 slice
  } else {
    kernelDims.z = kernelMatDims[2]; // slice = z-coordinate
  }

  // allocate memory for output
  plhs[0] = (mxArray *)mxCreateNumericArray(imageNDims, 
					    imageMatDims,
					    mxSINGLE_CLASS,
					    mxREAL);
  if (plhs[0] == NULL) {
    mexErrMsgTxt("Cannot allocate memory for output image");
  }

  // pointer to the Matlab input and output buffers
  float *convolvedImage =  (float *)mxGetData(prhs[1]);
  float *kernelImage =  (float *)mxGetData(prhs[2]);
  float *deconvolvedImage =  (float *)mxGetData(plhs[0]);
  
  // initialize Clarity by registering this application as a client
  Clarity_Register();

  // defaults for input parameters
  int iterations = 10;
  float epsilon = 0.01;

  // if input parameter is provided
  if (nrhs > 3) {

    // check that the parameter is given as a double
    if (!mxIsDouble(prhs[3])) {
      mexErrMsgTxt("Deconvolution parameter must be given as a double");
    }
    double *param = mxGetPr(prhs[3]);
    if (param == NULL) {
      mexErrMsgTxt("Couldn't read PARAM argument");
    }

    // epsilon parameter
    if (!strcmp(algorithmName, "w")) { // Wiener

      epsilon = (float)param[0];

    // iterations parameter
    } else if ((!strcmp(algorithmName, "jvc")) // JansenVanCittert
	       || (!strcmp(algorithmName, "ml"))) { // MaximumLikelihood

      iterations = (int)param[0];

    }
  }


  // run the deconvolution
  if (!strcmp(algorithmName, "w")) { // Wiener
    Clarity_WienerDeconvolve(convolvedImage, imageDims,
			     kernelImage, kernelDims, 
			     deconvolvedImage, epsilon);
  } else if (!strcmp(algorithmName, "jvc")) { // JansenVanCittert
    Clarity_JansenVanCittertDeconvolve(convolvedImage, imageDims,
				       kernelImage, kernelDims, 
				       deconvolvedImage, iterations);
  } else if (!strcmp(algorithmName, "ml")) { // MaximumLikelihood
    Clarity_MaximumLikelihoodDeconvolve(convolvedImage, imageDims,
					kernelImage, kernelDims, 
					deconvolvedImage, iterations);
  } else {
    mexErrMsgTxt("Deconvolution algorithm not implemented");
  }

  // unregister this application as a client
  Clarity_UnRegister();
  
  // successful return
  return;

}
