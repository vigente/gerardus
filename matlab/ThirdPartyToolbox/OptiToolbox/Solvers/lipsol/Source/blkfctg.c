/* C-Wrapper for LIPSOL FORTRAN Routines
 * J.Currie AUT May 2013
 */

#include "mex.h"

extern void BLKFCT(int *neqns, int *nsuper, int *xsuper, int *snode, int *split,
                   int *xlindx, int *lindx, int *xlnz, double *lnz, int *iwsiz,
                   int *iwork, int *tmpsiz, double *tmpvec, int *iflag, 
                   void (*MMPYN)(int*,int*,int*,int*,double*,double*,int*),
                   void (*SMXPY)(int*,int*,double*,int*,double*));

extern void MMPY1(int *m, int *n, int *q, int *xpnt, double *x, double *y, int *ldy);
extern void MMPY2(int *m, int *n, int *q, int *xpnt, double *x, double *y, int *ldy);
extern void MMPY4(int *m, int *n, int *q, int *xpnt, double *x, double *y, int *ldy);
extern void MMPY8(int *m, int *n, int *q, int *xpnt, double *x, double *y, int *ldy);
extern void SMXPY1(int *m, int *n, double *y, int *apnt, double *a);
extern void SMXPY2(int *m, int *n, double *y, int *apnt, double *a);
extern void SMXPY4(int *m, int *n, double *y, int *apnt, double *a);
extern void SMXPY8(int *m, int *n, double *y, int *apnt, double *a);
        
void mexFunction( int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[] )
{    
    int neqns, maxsub, nsuper, tmpsiz, level, flag, iwsiz, i;
    int *Pxlnz, *Pxsuper, *Psnode, *Psplit, *Pxlindx, *Plindx, *Piwork;
    double *PMxlnz, *PMxsuper, *PMsnode, *PMsplit, *PMxlindx, *PMlindx, *Plnz;
    double *Ptmpvec;

	//Check Inputs
	if(nrhs < 4)
		mexErrMsgTxt("BLKFCT requires 9 input arguments");
	//Get Sizes
	nsuper = (int)(max(mxGetM(prhs[1]),mxGetN(prhs[1]))-1);
    neqns = (int)max(mxGetM(prhs[2]),mxGetN(prhs[2]));
    maxsub = (int)max(mxGetM(prhs[5]),mxGetN(prhs[5]));
    tmpsiz = (int)mxGetScalar(prhs[7]);
    level = (int)mxGetScalar(prhs[8]);
    iwsiz = 2*neqns + 2*nsuper;

	//Get Input Args
    PMxlnz    = mxGetPr(prhs[0]);
    PMxsuper  = mxGetPr(prhs[1]);
    PMsnode   = mxGetPr(prhs[2]);
    PMsplit   = mxGetPr(prhs[3]);
    PMxlindx  = mxGetPr(prhs[4]);
    PMlindx   = mxGetPr(prhs[5]);    	
    
    //Create Output Args
    plhs[0] = mxDuplicateArray(prhs[6]); //BLKFCT edits in place
    Plnz    = mxGetPr(plhs[0]);
	
	//Create Internal Memory
    Pxlnz   = (int*)mxCalloc(neqns+1,sizeof(int));
    Pxsuper = (int*)mxCalloc(nsuper+1, sizeof(int));
    Psnode  = (int*)mxCalloc(neqns,  sizeof(int));
    Psplit  = (int*)mxCalloc(neqns,  sizeof(int));
    Pxlindx = (int*)mxCalloc(nsuper+1,sizeof(int));
    Plindx  = (int*)mxCalloc(maxsub,sizeof(int));
    Piwork  = (int*)mxCalloc(iwsiz,sizeof(int));
    Ptmpvec = (double*)mxCalloc(tmpsiz,sizeof(double));
	
	//Convert Input Args to Integers
    for(i=0;i<neqns;i++){
        Pxlnz[i]  = (int)PMxlnz[i];
        Psnode[i] = (int)PMsnode[i];
        Psplit[i] = (int)PMsplit[i];
    }
    Pxlnz[neqns] = (int)PMxlnz[neqns];
    for(i=0;i<maxsub;i++)
        Plindx[i] = (int)PMlindx[i];
    for(i=0;i<nsuper+1;i++) {
        Pxsuper[i]  = (int)PMxsuper[i];
        Pxlindx[i] = (int)PMxlindx[i];
    }    
	
    //Call FORTRAN BLKFCT (based on loop unrolling level)
    switch(level)
    {
        case 1:
            BLKFCT(&neqns,&nsuper,Pxsuper,Psnode,Psplit,Pxlindx,Plindx,Pxlnz,Plnz,&iwsiz,Piwork,&tmpsiz,Ptmpvec,&flag,MMPY1,SMXPY1);
            break;
        case 2:
            BLKFCT(&neqns,&nsuper,Pxsuper,Psnode,Psplit,Pxlindx,Plindx,Pxlnz,Plnz,&iwsiz,Piwork,&tmpsiz,Ptmpvec,&flag,MMPY2,SMXPY2);
            break;
        case 4:
            BLKFCT(&neqns,&nsuper,Pxsuper,Psnode,Psplit,Pxlindx,Plindx,Pxlnz,Plnz,&iwsiz,Piwork,&tmpsiz,Ptmpvec,&flag,MMPY4,SMXPY4);
            break;
        default:
            BLKFCT(&neqns,&nsuper,Pxsuper,Psnode,Psplit,Pxlindx,Plindx,Pxlnz,Plnz,&iwsiz,Piwork,&tmpsiz,Ptmpvec,&flag,MMPY8,SMXPY8);
            break;
    }

    //Check Init Status
    if(flag==-1)
        mexErrMsgTxt("Insufficient working storage in blkfct");
    
	//Free Internal Memory
	mxFree(Pxlnz);
	mxFree(Pxsuper);
	mxFree(Psnode);
	mxFree(Psplit);
    mxFree(Pxlindx);
    mxFree(Plindx);
    mxFree(Piwork);
    mxFree(Ptmpvec);
}