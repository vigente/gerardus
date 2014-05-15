/* QSOPTMEX - A MATLAB MEX Interface to QSOPT
 * Released Under the BSD 3-Clause License:
 * http://www.i2c2.aut.ac.nz/Wiki/OPTI/index.php/DL/License
 *
 * Copyright (C) Jonathan Currie 2013
 * www.i2c2.aut.ac.nz
 */

#include "mex.h"
#include "qsopt.h"
#include "reporter.h"
#include <exception>

#define QSOPT_VERSION "080725"

using namespace std;

//Function Prototypes
void printSolverInfo();
void checkInputs(const mxArray *prhs[], int nrhs);
double getStatus(int stat);

//Mex Printer (Reporter in QSOPT Terminology)
int mexPrinter(void *dest, const char *s)
{
    if(s != NULL) {
    	mexPrintf("%s",s);
        mexEvalString("drawnow;"); //flush draw buffer
    }
    return 0;   
}

void mexFunction(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[])
{
    //Input Args
    double *f, *A, *b, *lb, *ub;
    
    //Return Args
    double *x, *fval, *exitflag, *pi, *rc;
    
    //Options
    int maxiter, printLevel;
    double maxtime;
    
    //Internal Vars
    double *llb, *lub;
    size_t ncon, nin, ndec;
    size_t i;
    int rval, status, no = 0, a_lb = 0, a_ub = 0;
    mwIndex nzA;
    int *colCnt, *colBeg, *rowInd;
    char *sense;
    
    //Sparse Indicing
    mwIndex *A_ir, *A_jc;
    int *rows = NULL;
    //CoinBigIndex *cols = NULL;
    
    if(nrhs < 1) {
        if(nlhs < 1)
            printSolverInfo();
        else
            plhs[0] = mxCreateString(QSOPT_VERSION);
        return;
    }        
    
    //Check Inputs
    checkInputs(prhs,nrhs); 
    
    //Get pointers to Input variables
	f = mxGetPr(prhs[0]);
	A = mxGetPr(prhs[1]); 
    A_ir = mxGetIr(prhs[1]);
    A_jc = mxGetJc(prhs[1]);
    b = mxGetPr(prhs[2]);
    lb = mxGetPr(prhs[3]); 
    ub = mxGetPr(prhs[4]);

    //Get sizes
    ndec = mxGetM(prhs[0]);
    ncon = mxGetM(prhs[1]); 
    nin = (int)*mxGetPr(mxGetField(prhs[5],0,"nin"));
    nzA = A_jc[ndec];
    
    //Get options
    maxiter = (int)*mxGetPr(mxGetField(prhs[5],0,"maxiter"));
    printLevel = (int)*mxGetPr(mxGetField(prhs[5],0,"display"));
    maxtime = *mxGetPr(mxGetField(prhs[5],0,"maxtime"));
    
    //Create Outputs
    plhs[0] = mxCreateDoubleMatrix(ndec,1, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(1,1, mxREAL);
    plhs[2] = mxCreateDoubleMatrix(1,1, mxREAL);
    plhs[3] = mxCreateDoubleMatrix(ncon,1, mxREAL);
    plhs[4] = mxCreateDoubleMatrix(ndec,1, mxREAL);
    x = mxGetPr(plhs[0]); 
    fval = mxGetPr(plhs[1]); 
    exitflag = mxGetPr(plhs[2]);    
    pi = mxGetPr(plhs[3]);
    rc = mxGetPr(plhs[4]);
    

    //QSOPT Objects
    QSprob p;      

    //Allocate Vectors    
    colCnt = (int*)mxCalloc(ndec,sizeof(int));
    colBeg = (int*)mxCalloc(ndec,sizeof(int));
    rowInd = (int*)mxCalloc(nzA,sizeof(int));
    sense = (char*)mxCalloc(ncon,sizeof(char));

    //Process Bounds (must copy in as we will change them!)
    llb = (double*)mxCalloc(ndec,sizeof(double));
    lub = (double*)mxCalloc(ndec,sizeof(double));   
    //Create bounds if empty
    if(mxIsEmpty(prhs[3])) {
        for(i=0;i<ndec;i++)
            llb[i] = -QS_MAXDOUBLE;
    }
    else
        memcpy(llb,lb,ndec*sizeof(double));
    
    if(mxIsEmpty(prhs[4])) {
        for(i=0;i<ndec;i++)
            lub[i] = QS_MAXDOUBLE;
    }
    else
        memcpy(lub,ub,ndec*sizeof(double));
        
    //Ensure 'finite' bounds
    for(i = 0; i < ndec; i++) {
        if(mxIsInf(llb[i]))
            llb[i] = -QS_MAXDOUBLE;
        if(mxIsInf(lub[i]))
            lub[i] = QS_MAXDOUBLE;        
    }

    //Add Linear Constraints (Sparse!)
    if(ncon) {
        //Build Col Arrays
        for(i = 0; i < ndec; i++) {
            colCnt[i] = (int)(A_jc[i+1]-A_jc[i]);
            colBeg[i] = (int)A_jc[i];
        }
        //Build Row Array
        for(i = 0; i < nzA; i++)
            rowInd[i] = (int)A_ir[i];
        //Build Sense Array
        for(i = 0; i < ncon; i++) {
            if(i < nin)
                sense[i] = 'L';
            else
                sense[i] = 'E';
        }
    }
    else
        mexErrMsgTxt("You must supply A,b constraints to QSopt"); //just bounds

    //Load Problem into solver
    p = QSload_prob ("opti", (int)ndec, (int)ncon, colCnt, colBeg, rowInd, A, QS_MIN, f, b, sense, llb, lub, NULL, NULL);
    if(p != NULL)     
    {
        //Set Options        
        QSset_param(p,QS_PARAM_SIMPLEX_MAX_ITERATIONS, (int)maxiter);
        QSset_param_double(p,QS_PARAM_SIMPLEX_MAX_TIME, (double)maxtime); //doesn't seem to work
        
        //Setup Printing
        if(printLevel > 0) {
            QSset_param(p,QS_PARAM_SIMPLEX_DISPLAY, 1); //enable printing
            QSset_reporter(p, 10, mexPrinter, NULL); //iter skip seems to have a mind of its own...
        }
        
        //Solve Problem
        rval = QSopt_dual (p, &status);
        if(rval) { //failed
            *exitflag = -2;
            return;
        }
        //Solver ok - get solution
        QSget_solution (p, fval, x, pi, NULL, rc);  
        //Get Status
        *exitflag = getStatus(status);
    }
    else
        mexErrMsgTxt("Error Loading QSopt Problem!");

    //Clean up memory
    mxFree(colCnt);
    mxFree(colBeg);
    mxFree(rowInd);
    mxFree(sense);   
    mxFree(llb);
    mxFree(lub);
}               


//Check all inputs for size and type errors
void checkInputs(const mxArray *prhs[], int nrhs)
{
    size_t ndec, ncon;
    
    //Correct number of inputs
    if(nrhs < 6)
        mexErrMsgTxt("You must supply at least 6 arguments to qsopt (f, A, b, lb, ub, opts)"); 
    
    //Check we have an objective
    if(mxIsEmpty(prhs[0]))
        mexErrMsgTxt("You must supply an objective function!");
    
    //Check we have some constraints
    if(mxIsEmpty(prhs[1]) && mxIsEmpty(prhs[3]) && mxIsEmpty(prhs[5]) && mxIsEmpty(prhs[6]))
        mexErrMsgTxt("You have not supplied any constraints!");
   
    //Check options is a structure with correct fields
    if(!mxIsStruct(prhs[5]))
        mexErrMsgTxt("The options argument must be a structure!");
    else {
        if(mxGetFieldNumber(prhs[5],"maxiter") < 0)
            mexErrMsgTxt("The options structure should contain the field 'maxiter'");
        if(mxGetFieldNumber(prhs[5],"maxtime") < 0)
            mexErrMsgTxt("The options structure should contain the field 'maxtime'");
        if(mxGetFieldNumber(prhs[5],"nin") < 0)
            mexErrMsgTxt("The options structure should contain the field 'nin'");
        if(mxGetFieldNumber(prhs[5],"display") < 0)
            mexErrMsgTxt("The options structure should contain the field 'display'");
    }
    
    //Get Sizes
    ndec = mxGetM(prhs[0]);
    ncon = mxGetM(prhs[1]);
    
    //Check Constraint Pairs
    if(ncon && mxIsEmpty(prhs[2]))
        mexErrMsgTxt("b is empty!");
    
    //Check Sparsity (only supported in A)
    if(mxIsSparse(prhs[0]))
        mexErrMsgTxt("Only A is a sparse matrix");
    if(!mxIsSparse(prhs[1]))
        mexErrMsgTxt("A must be a sparse matrix");
    
    //Check Orientation
    if(mxGetM(prhs[0]) < mxGetN(prhs[0]))
        mexErrMsgTxt("f must be a column vector");
    if(mxGetM(prhs[2]) < mxGetN(prhs[2]))
        mexErrMsgTxt("b must be a column vector");
    if(mxGetM(prhs[3]) < mxGetN(prhs[3]))
        mexErrMsgTxt("lb must be a column vector");
    if(mxGetM(prhs[4]) < mxGetN(prhs[4]))
        mexErrMsgTxt("ub must be a column vector");
    
    //Check Sizes
    if(ncon) {
        if(mxGetN(prhs[1]) != ndec)
            mexErrMsgTxt("A has incompatible dimensions");
        if(mxGetM(prhs[2]) != ncon)
            mexErrMsgTxt("b has incompatible dimensions");
    }
    if(!mxIsEmpty(prhs[3]) && (mxGetM(prhs[3]) != ndec))
        mexErrMsgTxt("lb has incompatible dimensions");
    if(!mxIsEmpty(prhs[4]) && (mxGetM(prhs[4]) != ndec))
        mexErrMsgTxt("ub has incompatible dimensions");    
}


//Convert exiflag to OPTI status
double getStatus(int stat)
{
    double ret = -1;
    
    switch(stat)
    {
        case 1: //optimal
            ret = 1;
            break;
        case 2: //infeasible
            ret = -1;
            break;
        case 3: //max iterations
            ret = 0;
            break;
        case 4: //max time
            ret = 0;
            break;
        default: //inaccuracy / unbounded
            ret = -2;
            break;
    }
    return ret;
}  

//Print Solver Information
void printSolverInfo()
{    
    mexPrintf("\n-----------------------------------------------------------\n");
    mexPrintf(" QSOPT: Linear Programming Solver [v%s]\n",QSOPT_VERSION);      
    mexPrintf("  - Released under an Academic Only License: http://www2.isye.gatech.edu/~wcook/qsopt\n");
    mexPrintf("  - This is a closed source project\n");
    
    mexPrintf("\n MEX Interface J.Currie 2013 [BSD3] (www.i2c2.aut.ac.nz)\n");
    mexPrintf("-----------------------------------------------------------\n");
}
