/*
 *
 *
 * FORWARD_TV_mex Total variation of a 3D image
 *   This function should only be called by forward_TV.m
 
 
 * Author: Darryl McClymont <darryl.mcclymont@gmail.com>
 * Copyright © 2014 University of Oxford
 * Version: 0.1.1
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

#include <cmath> 
//#include <vector> 
#include <mex.h> 
#include <math.h> 
using namespace std; 


  
// Main function 
void mexFunction( int nlhs, mxArray *plhs[], 
                  int nrhs, const mxArray *prhs[] ) 
{ 
    double *INPUT, *DX, *DY, *DZ, *TV; 
    mwSize     ndim, ndim_zero; 
    const mwSize *dims, *dims_zero; 
  
  
    /* The input must be a noncomplex scalar double.*/
    ndim = mxGetNumberOfDimensions(prhs[0]);   
    dims = mxGetDimensions(prhs[0]); 
  
    
    // I couldn't work out how to get a single value output, so I'm just using the second input argument!
    ndim_zero = mxGetNumberOfDimensions(prhs[1]);   
    dims_zero = mxGetDimensions(prhs[1]); 
    
    /* Create matrices for the return arguments. */
    plhs[0] = mxCreateNumericArray(ndim, dims, mxDOUBLE_CLASS, mxREAL); 
    plhs[1] = mxCreateNumericArray(ndim, dims, mxDOUBLE_CLASS, mxREAL); 
    plhs[2] = mxCreateNumericArray(ndim, dims, mxDOUBLE_CLASS, mxREAL); 
    plhs[3] = mxCreateNumericArray(ndim_zero, dims_zero, mxDOUBLE_CLASS, mxREAL); 
  
    /* Assign pointers to each input and output. */
    DX = mxGetPr(plhs[0]); 
    DY = mxGetPr(plhs[1]); 
    DZ = mxGetPr(plhs[2]); 
    TV = mxGetPr(plhs[3]);
    INPUT = mxGetPr(prhs[0]); 
  
    
    int dim0 = (int)dims[0];
    int dim1 = (int)dims[1];
    int dim2 = (int)dims[2];
    

    
    int numel = dim0*dim1*dim2;
    
    
    // I had an implementation with lots of calls to sub2ind, but it was far too slow
    
    
    // Finite differences in the row dimension
    for (int i = 0; i < numel - 1; i++) 
    {  
        DX[i] = INPUT[i+1] - INPUT[i];
    }

    // Column dimension
    for (int i = 0; i < numel - dim0; i++) 
    {  
        DY[i] = INPUT[i + dim0] - INPUT[i];
    }
    
    // Slice dimension
    for (int i = 0; i < numel - dim0*dim1; i++) 
    {  
        DZ[i] = INPUT[i + dim0*dim1] - INPUT[i];
    }
    
    // Due to the indexing, the last row and column have values, but they should be zeros
    // DZ is fine though - the last slice is already zeros

    // fill the last row with zeros
    for (int i = dim0-1; i < numel; i+= dim0) 
    {  
        DX[i] = 0;
    }
    
    // do the same with the last column - this one is a bit trickier with the indexing
    int base_i = dim0 * (dim1-1);
    int k_index = 0;
    for (int k = 0; k < dim2; k++) 
    {  
        k_index = k * dim0 * dim1;
        for (int i = 0; i < dims[0]; i++)
        {
            // k index gets the slice, base_i gets us to the first row of the final column, and i gets us the rest of the way
            DY[k_index + base_i + i] = 0;
        }
    }
   
    
     // add them up to get the total variation
    double TV_adder = 0;
    
    for (int i = 0; i < numel; i++) 
    {  
       TV_adder+= abs(DX[i])+abs(DY[i])+abs(DZ[i]);
    }
        
    TV[0] = TV_adder;
}