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
C Usage: [perm, invp] = ordmmd(P) where diag(P) = 0
C
      SUBROUTINE MEXFUNCTION(NLHS, PLHS, NRHS, PRHS)
      INTEGER*4 PLHS(*), PRHS(*)
      INTEGER NLHS, NRHS
C
      INTEGER*4 MXCREATEFULL, MXGETPR
      INTEGER*4 MXCALLOC, MXGETIR, MXGETJC
      INTEGER MXGETM, MXGETN, MXGETNZMAX
C
C KEEP THE ABOVE SUBROUTINE, ARGUMENT, AND FUNCTION DECLARATIONS 
C FOR USE IN ALL YOUR FORTRAN MEX FILES.
C---------------------------------------------------------------------
C
      INTEGER*4 PMxadj, PMadjncy, PMperm, PMinvp
      INTEGER*4  Pxadj,  Padjncy,  Pperm,  Pinvp
      INTEGER*4 Piwork
      INTEGER neqns, iwsiz, flag, nofsub
C
C CHECK FOR PROPER NUMBER OF ARGUMENTS
C
      IF (NRHS .NE. 1) THEN
        CALL MEXERRMSGTXT('ORDMMD requires 1 input arguments')
      ELSEIF (NLHS .NE. 2) THEN
        CALL MEXERRMSGTXT('ORDMMD requires 2 output arguments')
      ENDIF
C
C CHECK FOR DIMENSIONS OF INPUT ARGUMENTS
C
      IF (MXISSPARSE(PRHS(1)) .NE. 1) THEN
        CALL MEXERRMSGTXT('Input matrix must be sparse')
      ENDIF
      neqns = MXGETM(PRHS(1))
      IF (neqns .NE. MXGETN(PRHS(1))) THEN
        CALL MEXERRMSGTXT('Input matrix must be square')
      ENDIF
C
C DEREFERENCE INPUT ARGUMENTS TO GET ARRAY POINTERS
C
      PMxadj   = MXGETJC(PRHS(1))
      PMadjncy = MXGETIR(PRHS(1))
C
C SPECIFY THE DIMENSIION OF WORKING VECTORS
C
      nzmax  = MXGETNZMAX(PRHS(1))
      iwsiz  = 4*neqns
C
C CREATE WORKING PARAMETERS
C
      Pxadj   = MXCALLOC(neqns+1,4)
      Padjncy = MXCALLOC(nzmax,  4)
      Pperm   = MXCALLOC(neqns,  4)
      Pinvp   = MXCALLOC(neqns,  4)
      Piwork  = MXCALLOC(iwsiz,  4)
C
C INPUT DATA CONVERSION
C
      call convertin(neqns, nzmax,
     & %VAL(PMxadj), %VAL(PMadjncy), 
     & %VAL( Pxadj), %VAL( Padjncy))
C
C DO THE ACTUAL COMPUTATIONS IN A FORTRAN SUBROUTINE
C
      call ordmmd(neqns, 
     & %VAL(Pxadj), %VAL(Padjncy), %VAL(Pinvp), %VAL(Pperm),
     & iwsiz, %VAL(Piwork), nofsub, flag)
C
C CHECK ERROR FLAG
C
      if (flag .eq. -1) then
        CALL MEXERRMSGTXT('Insufficient working storage 
     &                     in ordmmd.')
      endif
C
C CREATE MATRICES FOR RETURN ARGUMENTS
C
      PLHS(1) = MXCREATEDOUBLEMATRIX(neqns,1,0)
      PLHS(2) = MXCREATEDOUBLEMATRIX(neqns,1,0)
C
C DEREFERENCE OUTPUT ARGUMENTS TO GET REAL PART POINTERS
C
      PMperm   = MXGETPR(PLHS(1))
      PMinvp   = MXGETPR(PLHS(2))

      call int2real(neqns, 
     & %VAL( Pperm),%VAL( Pinvp),
     & %VAL(PMperm),%VAL(PMinvp))

      RETURN
      END


c------------------------------------------------------
c Convert mxadj and madjncy of range [0:n-1]
c        to  xadj and  adjncy of range [1:n];
c------------------------------------------------------
      subroutine convertin(neqns, nzmax, 
     &           mxadj, madjncy, 
     &            xadj,  adjncy )
      integer   neqns, nzmax, i
      integer   mxadj(neqns+1), madjncy(nzmax),
     &           xadj(neqns+1),  adjncy(nzmax)

      do 10 i = 1, neqns+1
         xadj(i) = mxadj(i) + 1
 10   continue

      do 20 i = 1, nzmax
         adjncy(i) = madjncy(i) + 1
 20   continue
      return
      end

c-------------------------------------------------
c Convert outputs from integer type in Fortran
c                   to real*8  type in Matlab
c-------------------------------------------------
      subroutine int2real(neqns,
     &            perm,  invp,
     &           mperm, minvp)

      integer  neqns
      integer   perm(1),  invp(1)
      real*8   mperm(1), minvp(1)

      do 10 i = 1, neqns
         mperm(i)  = perm(i)
         minvp(i)  = invp(i)
 10   continue
      return
      end
