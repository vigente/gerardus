/* C-Wrapper for LIPSOL FORTRAN Routines
 * J.Currie AUT May 2013
 */

#include "mex.h"

extern void INPNV(int *neqns, int *xadjf, int *adjf, double *anzf, int *perm, int *invp,
                  int *nsuper, int *xsuper, int *xlindx, int *lindx, int *xlnz,
                  double *lnz, int *offset);
        
void mexFunction( int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[] )
{    
    mwIndex *PMip, *PMjp;
    int neqns, nsuper, maxsub, nnzl, nzmax, i, j;
    int *Pip, *Pjp, *Pinvp, *Pperm, *Pxlnz, *Pxsuper, *Pxlindx, *Plindx, *Pxadjf, *Plink, *Padjf;
	double *Panzf, *PPdiag, *PP, *PMinvp, *PMperm, *PMxlnz, *PMxsuper, *PMxlindx, *PMlindx, *Plnz; 

	//Check Inputs
	if(nrhs < 9)
		mexErrMsgTxt("INPNV requires 9 input arguments");
	//Get Sizes
	neqns  = (int)max(mxGetM(prhs[1]),mxGetN(prhs[1]));
    nsuper = (int)(max(mxGetM(prhs[5]),mxGetN(prhs[5]))-1);
    maxsub = (int)max(mxGetM(prhs[7]),mxGetN(prhs[7]));
    nnzl   = (int)mxGetScalar(prhs[8]);
    nzmax  = (int)mxGetNzmax(prhs[1]);
	//Get Input Args
    PPdiag   = mxGetPr(prhs[0]);
    PP       = mxGetPr(prhs[1]);
    PMip     = mxGetIr(prhs[1]);
    PMjp     = mxGetJc(prhs[1]);
    PMinvp   = mxGetPr(prhs[2]);
    PMperm   = mxGetPr(prhs[3]);
    PMxlnz   = mxGetPr(prhs[4]);
    PMxsuper = mxGetPr(prhs[5]);
    PMxlindx = mxGetPr(prhs[6]);
    PMlindx  = mxGetPr(prhs[7]);
		
	//Create Outputs
	Plnz = mxGetPr(plhs[0] = mxCreateDoubleMatrix(nnzl,1,mxREAL));
	
	//Create Internal Memory
	Pip     = (int*)mxCalloc(nzmax, sizeof(int));
	Pjp     = (int*)mxCalloc(neqns+1, sizeof(int));
	Pinvp   = (int*)mxCalloc(neqns, sizeof(int));
	Pperm   = (int*)mxCalloc(neqns, sizeof(int));
    Pxlnz   = (int*)mxCalloc(neqns+1, sizeof(int));
    Pxsuper = (int*)mxCalloc(nsuper+1, sizeof(int));
    Pxlindx = (int*)mxCalloc(nsuper+1, sizeof(int));
    Plindx  = (int*)mxCalloc(maxsub, sizeof(int));
    Pxadjf  = (int*)mxCalloc(neqns+1, sizeof(int));
    Plink   = (int*)mxCalloc(neqns, sizeof(int));
    Padjf   = (int*)mxCalloc(nzmax+neqns, sizeof(int));
    Panzf   = (double*)mxCalloc(nzmax+neqns, sizeof(double));
	
	//Convert Input Args to Integers and FORTRAN Indices
	for(i=0;i<neqns;i++){
		Pjp[i] = (int)PMjp[i] + 1;
        Pinvp[i] = (int)PMinvp[i];
        Pperm[i] = (int)PMperm[i];
        Pxlnz[i] = (int)PMxlnz[i];        
    }
    Pjp[neqns] = (int)PMjp[neqns] + 1;
    Pxlnz[neqns] = (int)PMxlnz[neqns];
	for(i=0;i<nsuper+1;i++) {
		Pxsuper[i] = (int)PMxsuper[i];
        Pxlindx[i] = (int)PMxlindx[i];
    }
    for(i=0;i<maxsub;i++)
        Plindx[i] = (int)PMlindx[i];
    for(i=0;i<nzmax;i++)
        Pip[i] = (int)PMip[i] + 1;
    
    //Addiag Routine
    for(i=0;i<neqns;i++) {
        Panzf[Pjp[i]+i-1] = PPdiag[i];
        Padjf[Pjp[i]+i-1] = i+1;
        for(j=Pjp[i];j<Pjp[i+1];j++) {
            Panzf[j+i] = PP[j-1];
            Padjf[j+i] = Pip[j-1];
        }
        Pxadjf[i] = Pjp[i]+i;
    }
    Pxadjf[neqns] = Pjp[neqns]+neqns;
	
	//Call FORTRAN INPNV
	INPNV(&neqns, Pxadjf, Padjf, Panzf, Pperm, Pinvp, &nsuper, Pxsuper, Pxlindx, Plindx, Pxlnz, Plnz, Plink);
	
	//Free Internal Memory
	mxFree(Pip);
    mxFree(Pjp);
    mxFree(Pinvp);
    mxFree(Pperm);
    mxFree(Pxlnz);
	mxFree(Pxsuper);
	mxFree(Pxlindx);
	mxFree(Plindx);
    mxFree(Pxadjf);
    mxFree(Plink);
    mxFree(Padjf);
    mxFree(Panzf);
}