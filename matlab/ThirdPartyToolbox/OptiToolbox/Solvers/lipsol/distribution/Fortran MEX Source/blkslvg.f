C***********************************************************************
C
C   Version:        Beta-2.2
C   Last modified:  January 13, 1995
C   Author:         Yin Zhang
C                   Department of Mathematics and Statistics
C                   University of Maryland Baltimore County
C
C***********************************************************************
C
C Usage: [sol] = blkslv(xlnz,xsuper,xlindx,lindx,lnz,rrhs)
C
      SUBROUTINE MEXFUNCTION(NLHS, PLHS, NRHS, PRHS)
      INTEGER*4 PLHS(*), PRHS(*)
      INTEGER NLHS, NRHS
C
      INTEGER*4 MXCALLOC, MXGETPR
      INTEGER MXGETM, MXGETN
C
C KEEP THE ABOVE SUBROUTINE, ARGUMENT, AND FUNCTION DECLARATIONS 
C FOR USE IN ALL YOUR FORTRAN MEX FILES.
C---------------------------------------------------------------------
C
      INTEGER*4 PMxlnz, PMxsuper, PMxlindx, PMlindx
      INTEGER*4  Pxlnz,  Pxsuper,  Pxlindx,  Plindx
      INTEGER*4 Plnz, Prrhs, Psol
      INTEGER neqns, maxsub, nsuper
C
C CHECK FOR PROPER NUMBER OF ARGUMENTS
C
      IF (NRHS .NE. 6) THEN
        CALL MEXERRMSGTXT('SUPSLV requires 6 input arguments')
      ELSEIF (NLHS .NE. 1) THEN
        CALL MEXERRMSGTXT('SUPSLV requires 1 output argument')
      ENDIF
C
C SPECIFY THE DIMENSIION OF WORKING VECTORS
C
      neqns  = max(MXGETM(PRHS(6)),MXGETN(PRHS(6)))
      nsuper = max(MXGETM(PRHS(2)),MXGETN(PRHS(2))) - 1
      maxsub = max(MXGETM(PRHS(4)),MXGETN(PRHS(4)))
C
C DEREFERENCE ARGUMENTS TO GET ARRAY POINTERS
C
      PMxlnz   = MXGETPR(PRHS(1))
      PMxsuper = MXGETPR(PRHS(2))
      PMxlindx = MXGETPR(PRHS(3))
      PMlindx  = MXGETPR(PRHS(4))
      Plnz     = MXGETPR(PRHS(5))
      Prrhs    = MXGETPR(PRHS(6))
C
C CREATE WORKING PARAMETERS
C
      Pxlnz    = MXCALLOC(neqns+1, 4)
      Pxsuper  = MXCALLOC(nsuper+1, 4)
      Pxlindx  = MXCALLOC(nsuper+1, 4)
      Plindx   = MXCALLOC(maxsub,4 )
C
C INPUT DATA TRANSFORMATION
C
      call convertin(neqns, maxsub, nsuper,
     & %VAL(PMxlnz), %VAL(PMxsuper), 
     & %VAL( Pxlnz), %VAL( Pxsuper), 
     & %VAL(PMxlindx), %VAL(PMlindx),
     & %VAL( Pxlindx), %VAL( Plindx))
C
C CREATE OUTPUT PARAMETER
C
      PLHS(1) = MXCREATEDOUBLEMATRIX(neqns, 1, 0)
      Psol    = MXGETPR(PLHS(1))
C
C DO THE ACTUAL COMPUTATIONS IN A FORTRAN SUBROUTINE
C
      call blkslv(nsuper,
     & %VAL(Pxsuper), %VAL(Pxlindx), %VAL(Plindx),
     & %VAL(Pxlnz),   %VAL(Plnz),    %VAL(Prrhs))
C
C Copy solution to output vector
C
      call copysol(neqns, %VAL(Prrhs), %VAL(Psol))

      RETURN
      END


c----------------------------------------------------------------
c Convert from type real*8 in Matlab to type integer in Fortran
c----------------------------------------------------------------
      subroutine convertin(neqns, maxsub, nsuper,
     &           mxlnz, mxsuper, 
     &            xlnz,  xsuper, 
     &           mxlindx,mlindx,
     &            xlindx, lindx)

      integer neqns, maxsub, nsuper, i
      real*8  mxlnz(neqns+1), mxsuper(nsuper+1)
      integer  xlnz(neqns+1),  xsuper(nsuper+1)
      real*8  mxlindx(nsuper+1),mlindx(maxsub)
      integer  xlindx(nsuper+1), lindx(maxsub)


      do 10 i = 1, neqns + 1
         xlnz(i)   = mxlnz(i)
 10   continue

      do 20 i = 1, nsuper + 1
         xsuper(i) = mxsuper(i)
         xlindx(i) = mxlindx(i)
 20   continue

      do 30 i = 1, maxsub
         lindx(i) = mlindx(i)
 30   continue
      return
      end


c----------------------------------------------------------------
c Copy solution to output vector
c----------------------------------------------------------------
      subroutine copysol(neqns, rhs, sol)
      real*8 rhs(neqns), sol(neqns)
      integer i

      do 10 i = 1, neqns
         sol(i) = rhs(i)
 10   continue
      return
      end
