/* C-Wrapper for LIPSOL FORTRAN Routines
 * J.Currie AUT May 2013
 */

#include "mex.h"

extern void SFINIT(int *neqns, int *nnza, int *xadj, int *adjncy, int *perm,
                   int *invp, int *colcnt, int *nnzl, int *nsub, int *nsuper,
                   int *snode, int *xsuper, int *iwsiz, int *iwork, int *iflag);

extern void SYMFCT(int *neqns, int *adjlen, int *xadj, int *adjncy, int *perm,
                   int *invp, int *colcnt, int *nsuper, int *xsuper, int *snode,
                   int *nofsub, int *xlindx, int *lindx, int *xlnz, int *iwsiz,
                   int *iwork, int *iflag);

extern void BFINIT(int *neqns, int *nsuper, int *xsuper, int *snode, int *xlindx,
                   int *lindx, int *cachsz, int *tmpsiz, int *split);
        
void mexFunction( int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[] )
{    
    int neqns, nsub, nsuper, flag, nnzl, tmpsiz, iwsiz, anzmax, i, cachsz;
    int *Pxadj, *Padjncy, *Pperm, *Pinvp, *Pxlnz, *Pxsuper, *Plindx, *Pxlindx, *Psnode, *Psplit, *Pcolcnt, *Piwork;
    mwIndex *PMxadj, *PMadjncy;
    double *PMperm, *PMinvp, *PMxlnz, *PMnnzl, *PMxsuper, *PMxlindx, *PMlindx, *PMsnode, *PMsplit, *PMtmpsiz;

	//Check Inputs
	if(nrhs < 4)
		mexErrMsgTxt("SYMFCT requires 4 input arguments");
    if(nlhs < 8)
		mexErrMsgTxt("SYMFCT requires 8 output arguments");
	//Get Sizes
	neqns = (int)mxGetM(prhs[0]);   
    iwsiz = 7*neqns + 3;
    //Check Input Sizes
    if(mxGetN(prhs[0]) != neqns || !mxIsSparse(prhs[0]))
        mexErrMsgTxt("Input matrix must be a sparse and square");
    if(mxGetNumberOfElements(prhs[1]) != neqns)
        mexErrMsgTxt("PERM must be a neqn x 1 vector");
    if(mxGetNumberOfElements(prhs[2]) != neqns)
        mexErrMsgTxt("INVP must be a neqn x 1 vector");
    if(mxGetNumberOfElements(prhs[3]) != 1)
        mexErrMsgTxt("CACHSZ must be a scalar");

	//Get Input Args
    PMxadj   = mxGetJc(prhs[0]);
    PMadjncy = mxGetIr(prhs[0]);
    PMperm   = mxGetPr(prhs[1]);
    PMinvp   = mxGetPr(prhs[2]);
    cachsz   = (int)mxGetScalar(prhs[3]);
    anzmax   = (int)PMxadj[mxGetN(prhs[0])];	
	
	//Create Internal Memory
    Pxadj   = (int*)mxCalloc(neqns+1,sizeof(int));
    Padjncy = (int*)mxCalloc(anzmax, sizeof(int));
    Pperm   = (int*)mxCalloc(neqns,  sizeof(int));
    Pinvp   = (int*)mxCalloc(neqns,  sizeof(int));
    Pxlnz   = (int*)mxCalloc(neqns+1,sizeof(int));
    Pxsuper = (int*)mxCalloc(neqns+1,sizeof(int));
    Pxlindx = (int*)mxCalloc(neqns+1,sizeof(int));
    Psnode  = (int*)mxCalloc(neqns,  sizeof(int));
    Psplit  = (int*)mxCalloc(neqns,  sizeof(int));
    Pcolcnt = (int*)mxCalloc(neqns,  sizeof(int));
    Piwork  = (int*)mxCalloc(iwsiz,  sizeof(int));
	
	//Convert Input Args to Integers & To Fortran Indices
    for(i=0;i<neqns;i++){
        Pxadj[i] = (int)PMxadj[i] + 1;
        Pperm[i] = (int)PMperm[i];
        Pinvp[i] = (int)PMinvp[i];
    }
    Pxadj[neqns] = (int)PMxadj[neqns] + 1;
    for(i=0;i<anzmax;i++)
        Padjncy[i] = (int)PMadjncy[i] + 1;
	
    //Call FORTRAN SFINIT
    SFINIT(&neqns, &anzmax, Pxadj, Padjncy, Pperm, Pinvp, Pcolcnt, &nnzl, &nsub, &nsuper,
           Psnode, Pxsuper, &iwsiz, Piwork, &flag);
    
    //Check Init Status
    if(flag==-1)
        mexErrMsgTxt("Insufficient working storage in sfinit");
    
    //Create Outputs (need nsuper and nsub)
	PMxlnz   = mxGetPr(plhs[0] = mxCreateDoubleMatrix(neqns+1,1,mxREAL));
    PMnnzl   = mxGetPr(plhs[1] = mxCreateDoubleMatrix(1,1,mxREAL));
    PMxsuper = mxGetPr(plhs[2] = mxCreateDoubleMatrix(nsuper+1,1,mxREAL));
    PMxlindx = mxGetPr(plhs[3] = mxCreateDoubleMatrix(nsuper+1,1,mxREAL));
    PMlindx  = mxGetPr(plhs[4] = mxCreateDoubleMatrix(nsub,1,mxREAL));
    PMsnode  = mxGetPr(plhs[5] = mxCreateDoubleMatrix(neqns,1,mxREAL));
    PMsplit  = mxGetPr(plhs[6] = mxCreateDoubleMatrix(neqns,1,mxREAL));
    PMtmpsiz = mxGetPr(plhs[7] = mxCreateDoubleMatrix(1,1,mxREAL));
    
    //Create Further Working Memory
    Plindx = (int*)mxCalloc(2*nsub, sizeof(int));
    
    //Call FORTRAN SYMFCT
    SYMFCT(&neqns, &anzmax, Pxadj, Padjncy, Pperm, Pinvp, Pcolcnt, &nsuper, Pxsuper, Psnode,
           &nsub, Pxlindx, Plindx, Pxlnz, &iwsiz, Piwork, &flag);

    //Check symfct status
    if(flag==-1)
        mexErrMsgTxt("Insufficient working storage in symfct");
    if(flag==-2)
        mexErrMsgTxt("Inconsistency in the input for symfct");
    
    //Call FORTRAN BFINIT
    BFINIT(&neqns, &nsuper, Pxsuper, Psnode, Pxlindx, Plindx, &cachsz, &tmpsiz, Psplit);

    //Return Integer Results as Doubles
    *PMnnzl = (double)nnzl;
    *PMtmpsiz = (double)tmpsiz;
	for(i=0;i<neqns;i++) {
		PMxlnz[i] = (double)Pxlnz[i];
		PMsnode[i] = (double)Psnode[i];
        PMsplit[i] = (double)Psplit[i];
        PMperm[i] = (double)Pperm[i];
        PMinvp[i] = (double)Pinvp[i];
	}
    PMxlnz[neqns] = (double)Pxlnz[neqns];
    for(i=0;i<nsuper+1;i++) {
        PMxsuper[i] = (double)Pxsuper[i];
		PMxlindx[i] = (double)Pxlindx[i];
    }
    for(i=0;i<nsub;i++)
        PMlindx[i] = (double)Plindx[i];
    
	//Free Internal Memory
	mxFree(Pxadj);
	mxFree(Padjncy);
	mxFree(Pperm);
	mxFree(Pinvp);
    mxFree(Pxlnz);
    mxFree(Pxsuper);
    mxFree(Pxlindx);
    mxFree(Plindx);
    mxFree(Psnode);
    mxFree(Psplit);
    mxFree(Pcolcnt);
    mxFree(Piwork);
}