C***********************************************************************
C
C   Version:        Beta-2.2
C   Last modified:  February 14, 1995
C   Author:         Yin Zhang
C                   Department of Mathematics and Statistics
C                   University of Maryland Baltimore County
C
C***********************************************************************
C
C Usage: [lnz] = blkfct(xlnz,xsuper,snode,split,xlindx,lindx,...
C                       lnz,tmpsiz,level)
C
      SUBROUTINE MEXFUNCTION(NLHS, PLHS, NRHS, PRHS)
      INTEGER*4 PLHS(*), PRHS(*)
      INTEGER NLHS, NRHS
C
      INTEGER*4 MXGETPR, MXCALLOC
      INTEGER MXGETM, MXGETN
      REAL*8  MXGETSCALAR
C
C KEEP THE ABOVE SUBROUTINE, ARGUMENT, AND FUNCTION DECLARATIONS 
C FOR USE IN ALL YOUR FORTRAN MEX FILES.
C---------------------------------------------------------------------
C
      INTEGER*4 PMxlnz, PMxsuper
      INTEGER*4  Pxlnz,  Pxsuper
      INTEGER*4 PMsnode, PMsplit, PMxlindx,PMlindx
      INTEGER*4  Psnode,  Psplit,  Pxlindx, Plindx
      INTEGER*4 Plnz, Piwork, Ptmpvec
      INTEGER neqns, maxsub, nsuper, tmpsiz, level, flag, iwsiz
      EXTERNAL  mmpy1,  mmpy2,  mmpy4,  mmpy8
      EXTERNAL smxpy1, smxpy2, smxpy4, smxpy8
C
C CHECK FOR PROPER NUMBER OF ARGUMENTS
C
      IF (NRHS .NE. 9) THEN
        CALL MEXERRMSGTXT('BLKFCT requires 9 input arguments')
      ELSEIF (NLHS .NE. 1) THEN
        CALL MEXERRMSGTXT('BLKFCT requires 1 output arguments')
      ENDIF
C
C SPECIFY THE DIMENSIION OF WORKING VECTORS
C
      nsuper = max(MXGETM(PRHS(2)),MXGETN(PRHS(2))) - 1
      neqns  = max(MXGETM(PRHS(3)),MXGETN(PRHS(3)))
      maxsub = max(MXGETM(PRHS(6)),MXGETN(PRHS(6)))
      tmpsiz = MXGETSCALAR(PRHS(8))
      level  = MXGETSCALAR(PRHS(9))
      iwsiz  = 2*neqns + 2*nsuper
C
C DEREFERENCE ARGUMENTS TO GET ARRAY POINTERS
C
      PMxlnz    = MXGETPR(PRHS(1))
      PMxsuper  = MXGETPR(PRHS(2))
      PMsnode   = MXGETPR(PRHS(3))
      PMsplit   = MXGETPR(PRHS(4))
      PMxlindx  = MXGETPR(PRHS(5))
      PMlindx   = MXGETPR(PRHS(6))
      Plnz      = MXGETPR(PRHS(7))
C
C CREATE WORKING PARAMETERS
C
      Pxlnz   = MXCALLOC(neqns+1, 4)
      Pxsuper = MXCALLOC(nsuper+1,4)
      Psnode  = MXCALLOC(neqns,   4)
      Psplit  = MXCALLOC(neqns,   4)
      Pxlindx = MXCALLOC(nsuper+1,4)
      Plindx  = MXCALLOC(maxsub,  4)
      Piwork  = MXCALLOC(iwsiz,   4)
      Ptmpvec = MXCALLOC(tmpsiz,  8)
C
C INPUT DATA TRANSFORMATION
C
      call convertin(neqns, maxsub, nsuper,
     & %VAL(PMxlnz),  %VAL(PMxsuper), %VAL(PMsnode),
     & %VAL(PMsplit), %VAL(PMxlindx), %VAL(PMlindx),
     & %VAL( Pxlnz),  %VAL( Pxsuper), %VAL( Psnode),
     & %VAL( Psplit), %VAL( Pxlindx), %VAL( Plindx))
C
C CREATE OUTPUT PARAMETERS
C
      PLHS(1) = PRHS(7)
C
C DO THE ACTUAL COMPUTATIONS IN A FORTRAN SUBROUTINE
C
      if (level .eq. 1) then
        call blkfct(neqns, nsuper, %VAL(Pxsuper), %VAL(Psnode), 
     &  %VAL(Psplit), %VAL(Pxlindx), %VAL(Plindx), %VAL(Pxlnz), 
     &  %VAL(Plnz), iwsiz, %VAL(Piwork), tmpsiz, %VAL(Ptmpvec),
     &  flag, mmpy1, smxpy1)
      endif
      if (level .eq. 2) then
        call blkfct(neqns, nsuper, %VAL(Pxsuper), %VAL(Psnode), 
     &  %VAL(Psplit), %VAL(Pxlindx), %VAL(Plindx), %VAL(Pxlnz), 
     &  %VAL(Plnz), iwsiz, %VAL(Piwork), tmpsiz, %VAL(Ptmpvec),
     &  flag, mmpy2, smxpy2)
      endif
      if (level .eq. 4) then
        call blkfct(neqns, nsuper, %VAL(Pxsuper), %VAL(Psnode), 
     &  %VAL(Psplit), %VAL(Pxlindx), %VAL(Plindx), %VAL(Pxlnz), 
     &  %VAL(Plnz), iwsiz, %VAL(Piwork), tmpsiz, %VAL(Ptmpvec),
     &  flag, mmpy4, smxpy4)
      endif
      if (level .eq. 8) then
        call blkfct(neqns, nsuper, %VAL(Pxsuper), %VAL(Psnode), 
     &  %VAL(Psplit), %VAL(Pxlindx), %VAL(Plindx), %VAL(Pxlnz), 
     &  %VAL(Plnz), iwsiz, %VAL(Piwork), tmpsiz, %VAL(Ptmpvec),
     &  flag, mmpy8, smxpy8)
      endif
C
C CHECK ERROR FLAG
C
      if (flag .lt. -1) then
        call MEXERRMSGTXT('Insufficient working storage 
     &                     in blkfct.')
      endif

C
      RETURN
      END


c----------------------------------------------------------------
c Convert from type real*8 in Matlab to type integer in Fortran
c----------------------------------------------------------------
      subroutine convertin(neqns, maxsub, nsuper,
     &           mxlnz, mxsuper, msnode, msplit, mxlindx, mlindx,
     &            xlnz,  xsuper,  snode,  split,  xlindx,  lindx)
      integer neqns, maxsub, nsuper, i
      real*8  mxlnz(1), mxsuper(1)
      integer  xlnz(1),  xsuper(1)
      real*8  msnode(1),msplit(1),mxlindx(1),mlindx(1)
      integer  snode(1), split(1), xlindx(1), lindx(1)
    
      do 10 i = 1, neqns
         xlnz(i)  = mxlnz(i)
         snode(i) = msnode(i)
         split(i) = msplit(i)
 10   continue
      xlnz(neqns+1) = mxlnz(neqns+1)

      do 20 i = 1, maxsub
         lindx(i) = mlindx(i)
 20   continue
      do 30 i = 1, nsuper + 1
         xsuper(i) = mxsuper(i)
         xlindx(i) = mxlindx(i)
 30   continue
      return
      end
