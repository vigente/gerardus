/* C-Wrapper for LIPSOL FORTRAN Routines
 * J.Currie AUT May 2013
 */

#include "mex.h"

extern void ORDMMD(int *neqns, int *xadj, int *adjncy, int *invp, int *perm,
                   int *iwsiz, int *iwork, int *nofsub, int *iflag);
        
void mexFunction( int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[] )
{    
	mwIndex *PMxadj, *PMadjncy;	
    int neqns, nnz, flag, iwsiz, nofsub, i;
	int *Pxadj, *Padjncy, *Pperm, *Pinvp, *Piwork;
	double *PMperm, *PMinvp;

	//Check Inputs
	if(nrhs < 1)
		mexErrMsgTxt("ORDMMD requires 1 input argument");
	if(!mxIsSparse(prhs[0]) || !mxIsDouble(prhs[0]))
		mexErrMsgTxt("P must be a sparse, double matrix");
	//Get Size
	neqns = (int)mxGetM(prhs[0]);
	iwsiz = 4*neqns;
	//Check P is Square
	if(neqns != mxGetN(prhs[0]))
		mexErrMsgTxt("P must be a square matrix!");
		
	//Create Outputs
	PMperm = mxGetPr(plhs[0] = mxCreateDoubleMatrix(neqns,1,mxREAL));
    PMinvp = mxGetPr(plhs[1] = mxCreateDoubleMatrix(neqns,1,mxREAL));
	//Get Sparse Indices			
	PMxadj = mxGetJc(prhs[0]);
	PMadjncy = mxGetIr(prhs[0]);
	nnz = (int)PMxadj[neqns];
	
	//Create Internal Memory
	Pxadj   = (int*)mxCalloc(neqns+1, sizeof(int));
	Padjncy = (int*)mxCalloc(nnz, sizeof(int));
	Pperm   = (int*)mxCalloc(neqns, sizeof(int));
	Pinvp   = (int*)mxCalloc(neqns, sizeof(int));
	Piwork  = (int*)mxCalloc(iwsiz, sizeof(int));
	
	//Copy Cast Sparse Indices & Convert to Fortran Indices
	for(i=0;i<=neqns;i++)
		Pxadj[i] = (int)PMxadj[i] + 1;
	for(i=0;i<nnz;i++)
		Padjncy[i] = (int)PMadjncy[i] + 1;
	
	//Call FORTRAN ORDMMD
	ORDMMD(&neqns, Pxadj, Padjncy, Pinvp, Pperm, &iwsiz, Piwork, &nofsub, &flag);
	
	//Check Error Flag
	if(flag == -1)
		mexErrMsgTxt("Insufficient working storage in ordmmd");
	
	//Return Integer Results as Doubles
	for(i=0;i<neqns;i++) {
		PMperm[i] = (double)Pperm[i];
		PMinvp[i] = (double)Pinvp[i];
	}
	
	//Free Internal Memory
	mxFree(Pxadj);
	mxFree(Padjncy);
	mxFree(Pperm);
	mxFree(Pinvp);
	mxFree(Piwork);
}