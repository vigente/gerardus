/* C-Wrapper for LIPSOL FORTRAN Routines
 * J.Currie AUT May 2013
 */

#include "mex.h"

extern void BLKSLV(int *nsuper, int *xsuper, int *xlindx, int *lindx, int *xlnz,
                   double *lnz, double *rhs);
        
void mexFunction( int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[] )
{    
    int neqns, nsuper, maxsub, i;
    int *Pxlnz, *Pxsuper, *Pxlindx, *Plindx;
	double *PMxlnz, *PMxsuper, *PMxlindx, *PMlindx, *Plnz, *Prrhs, *Psol;

	//Check Inputs
	if(nrhs < 6)
		mexErrMsgTxt("BLKSLV requires 6 input arguments");
	//Get Sizes
	neqns = (int)max(mxGetM(prhs[5]),mxGetN(prhs[5]));
    nsuper = (int)(max(mxGetM(prhs[1]),mxGetN(prhs[1]))-1);
    maxsub = (int)max(mxGetM(prhs[3]),mxGetN(prhs[3]));
	//Get Input Args
    PMxlnz   = mxGetPr(prhs[0]);
    PMxsuper = mxGetPr(prhs[1]);
    PMxlindx = mxGetPr(prhs[2]);
    PMlindx  = mxGetPr(prhs[3]);
    Plnz     = mxGetPr(prhs[4]);
    Prrhs    = mxGetPr(prhs[5]);
		
	//Create Outputs
	Psol = mxGetPr(plhs[0] = mxCreateDoubleMatrix(neqns,1,mxREAL));
	
	//Create Internal Memory
	Pxlnz   = (int*)mxCalloc(neqns+1, sizeof(int));
	Pxsuper = (int*)mxCalloc(nsuper+1, sizeof(int));
	Pxlindx = (int*)mxCalloc(nsuper+1, sizeof(int));
	Plindx  = (int*)mxCalloc(maxsub, sizeof(int));
	
	//Convert Input Args to Integers
	for(i=0;i<neqns+1;i++)
		Pxlnz[i] = (int)PMxlnz[i];
	for(i=0;i<nsuper+1;i++) {
		Pxsuper[i] = (int)PMxsuper[i];
        Pxlindx[i] = (int)PMxlindx[i];
    }
    for(i=0;i<maxsub;i++)
        Plindx[i] = (int)PMlindx[i];
	
	//Call FORTRAN BLKSLV
	BLKSLV(&nsuper, Pxsuper, Pxlindx, Plindx, Pxlnz, Plnz, Prrhs);
	
	//Copy Solution
    memcpy(Psol,Prrhs,neqns*sizeof(double));
	
	//Free Internal Memory
	mxFree(Pxlnz);
	mxFree(Pxsuper);
	mxFree(Pxlindx);
	mxFree(Plindx);
}