/* Fast 3D linear interpolation 
   Andriy Myronenko, Feb 2008, email: myron@csee.ogi.edu, homepage: http://www.bme.ogi.edu/~myron/
 
Output_image = mirt3D_mexinterp(Input_image, XI,YI,ZI) interpolates the 3D image 'Input_image' at
   the points with coordinates XI,YI,ZI. Input_image is assumed to  be defined at a regular grid 1:N, 1:M, 1:K, 
   where [M,N,K]=size(Input_images). Points outside the boundary return NaNs. 
   This is equivalent to Matlab's:
   Output_image = interp3(Input_image,XI,YI,ZI,'linear',NaN);
 
Output_images = mirt3D_mexinterp(Input_images, XI,YI,ZI). Input_images can be a stack of many 3D images (4D).
  The function interpolates each of the 3D images at X,Y,Z coordinates and return a stack of corresponding
  interpolated images. This is equivalent to Matlabs

  Input_images(:,:,:,1)=Input_image1;
  Input_images(:,:,:,2)=Input_image2;
  Input_images(:,:,:,3)=Input_image3;
  Input_images(:,:,:,4)=Input_image4;
 
 Output_images(:,:,:,1) = interp3(Input_image1,XI,YI,ZI,'linear',NaN);
  Output_images(:,:,:,2) = interp3(Input_image2,XI,YI,ZI,'linear',NaN);
  Output_images(:,:,:,3) = interp3(Input_image3,XI,YI,ZI,'linear',NaN);
  Output_images(:,:,:,4) = interp3(Input_image4,XI,YI,ZI,'linear',NaN);

 This is especially usefull to vector valued 3D images, RGB images, or to interpolate the whole 3D video at the same coordinates
 or to interpolate image and its gradients at the same time (in image registration).
 The speed gain is from the precomputation of nearest points for interpolation, which are the same for all images in a stack.

*/


#include <math.h>
#include "mex.h"


void mirt3D_mexinterp(
double* Z,
double* S,
double* T,
double* W,
double* F,
int	MN,
int nrows,
int ncols,
int npages,
int ndim
)
{
     
    int	n, in1, in2, in3, in4, in5, in6, in7, in8;
    double	t, s,  s1, w, w1, tmp, nan;
    double m1, m2, m3, m4, m5, m6, m7, m8;
    int ndx, nw, Zshift, i, nrowsncols, ft, fs, fw;

    nw = nrows*ncols;
    nrowsncols=nw*npages;
    nan=mxGetNaN();
    
    for (n=0; n < MN; n++) {
        
        t=T[n];
        s=S[n];
        w=W[n];
           
        ft=(int) floor(t);
        fs=(int) floor(s);
        fw=(int) floor(w);
        
        
        if (fs<1 || s>ncols || ft<1 || t>nrows || fw<1 || w>npages){
             /* Put nans if outside*/
            for (i = 0; i < ndim; i++) F[n+i*MN]=nan; }
        else  {
            
         
            ndx =  ft+(fs-1)*nrows+(fw-1)*nw;
            
            if (s==ncols){ s=s+1; ndx=ndx-nrows; }
            s=s-fs;
            if (t==nrows){  t=t+1; ndx=ndx-1; }
            t=t-ft;
            if (w==npages){  w=w+1; ndx=ndx-nw; }
            w=w-fw;
           
            in1=ndx-1;
            in2=ndx;
            // in3=ndx+nrows-1;
            in3=in1+nrows;
            // in4=ndx+nrows;
            in4=in3+1;
            
            // in5=ndx+nw-1;
            in5=in1+nw;
           // in6=ndx+nw;
            in6=in5+1;
           // in7=ndx+nrows+nw-1;
            in7=in5+nrows;
           // in8=ndx+nrows+nw;
            in8 = in7+1;
            
          ////////////
          s1=1-s;
          w1=1-w;
                    
          tmp=s1*w1;
          m2=t*tmp;
          m1=tmp-m2;
                                        
          tmp=s*w1;
          m4=t*tmp;
          m3=tmp-m4;
                                       
          tmp=s1*w;
          m6=t*tmp;
          m5=tmp-m6;
                                        
          tmp=s*w;
          m8=t*tmp;
          m7=tmp-m8;


             for (i = 0; i < ndim; i++){
                   
                    Zshift=i*nrowsncols;
                    F[n+i*MN]=Z[in1+Zshift]*m1+Z[in2+Zshift]*m2+Z[in3+Zshift]*m3+Z[in4+Zshift]*m4+Z[in5+Zshift]*m5+Z[in6+Zshift]*m6+
                            Z[in7+Zshift]*m7+Z[in8+Zshift]*m8;
            }
 
        } 
   
    } // cycle end
   
    return;
}

    
    
   
// ------- Main function definitions ----------
/* Input arguments */
#define IN_Z		prhs[0]
#define IN_S		prhs[1]
#define IN_T		prhs[2]
#define IN_W		prhs[3]

/* Output arguments */
#define OUT_F		plhs[0]

/* Gateway routine */
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
    double *Z, *S, *T, *W, *F;
    int  i, MN, nrows, ncols, npages, vol, ndim, newXndim, Xndim, *newdims;
    const int *dims, *Xdims;
    
    
    /* Check for input errors */
    if (nlhs>1)
    mexErrMsgTxt("Wrong number of output parameters, usage:  Output_images = mirt3D_mexinterp(Input_images, X, Y, Z)");
    if (nrhs!=4)
    mexErrMsgTxt("Wrong number of input parameters, usage:  Output_images = mirt3D_mexinterp(Input_images, X, Y, Z)");

    if (!mxIsDouble(IN_Z) || !mxIsDouble(IN_S) || !mxIsDouble(IN_T) || !mxIsDouble(IN_W))
    mexErrMsgTxt("mirt3D_mexinterp: Input arguments must be double.");
    
    if ((mxGetNumberOfDimensions(IN_S) != mxGetNumberOfDimensions(IN_T)) ||
      (mxGetNumberOfDimensions(IN_S) != mxGetNumberOfDimensions(IN_W)) ||
      (mxGetNumberOfElements(IN_S) != mxGetNumberOfElements(IN_T)) ||
      (mxGetNumberOfElements(IN_S) != mxGetNumberOfElements(IN_W))) mexErrMsgTxt("Input parameters X, Y, Z must have the same size");
    
    /* Get the sizes of each input argument */
    
    Xndim = mxGetNumberOfDimensions(IN_S);
    Xdims = mxGetDimensions(IN_S);
    
    
    ndim = mxGetNumberOfDimensions(IN_Z);
    dims = mxGetDimensions(IN_Z);
    newdims = (int*) calloc(ndim-1, sizeof(int));
       
    MN=1;
    for (i = 0; i < Xndim; i++) {MN =MN*Xdims[i];};  /*Total number of interpolations points in 1 image*/
       
    
    vol=1; newXndim=Xndim;
      if (ndim>3) {   /*Check if we have several images*/
        
         
        if ((Xndim==2) && (Xdims[1]==1))  {newXndim=newXndim-1; }  /*Check if interpolate along column vectors*/
        newdims = (int*) calloc(newXndim+1, sizeof(int));           /*Allocate space for the new number of dimensions for output*/
        for (i = 0; i < newXndim; i++) {newdims[i]=Xdims[i];};  /*Copy original dimenstions*/
        newdims[newXndim]=dims[3];                             /*Add the number of  images as a last dimenstion*/
        newXndim=newXndim+1;                                       /*Set the new number of dimenstions*/ 
        vol=dims[3];}
    else
    {
       newdims = (int*) calloc(newXndim, sizeof(int));
       for (i = 0; i < newXndim; i++) {newdims[i]=Xdims[i];};
        
    }
    
   
    /*Create the array to put the interpolated points*/
    OUT_F = mxCreateNumericArray(newXndim, newdims, mxDOUBLE_CLASS, mxREAL);
    
    /* Input image size */
    nrows = dims[0];
    ncols = dims[1];
    npages = dims[2];
    
  /* Assign pointers to the input arguments */
    Z = mxGetPr(IN_Z);
    S = mxGetPr(IN_S);
    T = mxGetPr(IN_T);
    W = mxGetPr(IN_W);
    
  /* Assign pointers to the output arguments */
    F      = mxGetPr(OUT_F);
    
  /* Do the actual computations in a subroutine */
    mirt3D_mexinterp(Z, S, T, W, F, MN, nrows, ncols, npages, vol);
    
    free((void*)newdims);
    return;
}


