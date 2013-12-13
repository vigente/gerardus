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
C Usage: [xlnz,nnzl,xsuper,xlindx,lindx,snode,split,tmpsiz,perm,invp]...
C                                     = symfct(P, perm, invp, cachsz)
C
      SUBROUTINE MEXFUNCTION(NLHS, PLHS, NRHS, PRHS)
      INTEGER*4 PLHS(*), PRHS(*)
      INTEGER NLHS, NRHS
C
      INTEGER*4 MXCREATEFULL, MXGETPR
      INTEGER*4 MXCALLOC, MXGETIR, MXGETJC
      INTEGER MXGETM, MXGETN, MXGETNZMAX
      REAL*8  MXGETSCALAR
C
C KEEP THE ABOVE SUBROUTINE, ARGUMENT, AND FUNCTION DECLARATIONS 
C FOR USE IN ALL YOUR FORTRAN MEX FILES.
C---------------------------------------------------------------------
C
      INTEGER*4 PMxadj, PMadjncy, PMperm, PMinvp, PMxlnz, PMxsuper
      INTEGER*4  Pxadj,  Padjncy,  Pperm,  Pinvp,  Pxlnz,  Pxsuper
      INTEGER*4 PMxlindx, PMlindx, PMsnode, PMsplit, PMnnzl, PMtmpsiz
      INTEGER*4  Pxlindx,  Plindx,  Psnode,  Psplit
      INTEGER*4 Pcolcnt, Piwork
      INTEGER m, n, cachsz, iwsiz, flag, nnzl, tmpsiz
      INTEGER neqns, anzmax, nsub, nsuper
C
C CHECK FOR PROPER NUMBER OF ARGUMENTS
C
      IF (NRHS .NE. 4) THEN
        CALL MEXERRMSGTXT('SYMFCT requires 4 input arguments')
      ELSEIF (NLHS .NE. 10) THEN
        CALL MEXERRMSGTXT('SYMFCT requires 10 output arguments')
      ENDIF
C
C CHECK FOR DIMENSIONS OF INPUT ARGUMENTS
C
      neqns = MXGETM(PRHS(1))
      IF (neqns .NE. MXGETN(PRHS(1))) THEN
        CALL MEXERRMSGTXT('Input matrix must be square')
      ENDIF
C
      m = MXGETM(PRHS(2))
      n = MXGETN(PRHS(2))
      IF ((MAX(m,n) .NE. neqns) .OR. (MIN(m,n) .NE. 1)) THEN
        CALL MEXERRMSGTXT('SYMFCT requires PERM to be neqns x 1')
      ENDIF
C
      m = MXGETM(PRHS(3))
      n = MXGETN(PRHS(3))
      IF ((MAX(m,n) .NE. neqns) .OR. (MIN(m,n) .NE. 1)) THEN
        CALL MEXERRMSGTXT('SYMFCT requires INVP to be neqns x 1')
      ENDIF
C
      m = MXGETM(PRHS(4))
      n = MXGETN(PRHS(4))
      IF ((m .NE. 1) .OR. (n .NE. 1)) THEN
        CALL MEXERRMSGTXT('SYMFCT requires CACHSZ to be 1 x 1')
      ENDIF
C
C DEREFERENCE INPUT ARGUMENTS TO GET ARRAY POINTERS
C
      PMxadj   = MXGETJC(PRHS(1))
      PMadjncy = MXGETIR(PRHS(1))
      PMperm   = MXGETPR(PRHS(2))
      PMinvp   = MXGETPR(PRHS(3))
C
C SPECIFY THE DIMENSIION OF WORKING VECTORS
C
      anzmax = MXGETNZMAX(PRHS(1))
      iwsiz  = 7*neqns + 3
C
C CREATE WORKING PARAMETERS
C
      Pxadj   = MXCALLOC(neqns+1,4)
      Padjncy = MXCALLOC(anzmax, 4)
      Pperm   = MXCALLOC(neqns,  4)
      Pinvp   = MXCALLOC(neqns,  4)
      Pxlnz   = MXCALLOC(neqns+1,4)
      Pxsuper = MXCALLOC(neqns+1,4)
      Pxlindx = MXCALLOC(neqns+1,4)
      Psnode  = MXCALLOC(neqns,  4)
      Psplit  = MXCALLOC(neqns,  4)
      Pcolcnt = MXCALLOC(neqns,  4)
      Piwork  = MXCALLOC(iwsiz,  4)
C
C INPUT DATA CONVERSION
C
      call convertin(neqns, anzmax,
     & %VAL(PMxadj), %VAL(PMadjncy), %VAL(PMperm), %VAL(PMinvp),
     & %VAL( Pxadj), %VAL( Padjncy), %VAL( Pperm), %VAL( Pinvp))
C
C CALL THE ACTUAL FORTRAN SUBROUTINES
C
      call sfinit(neqns, anzmax,
     & %VAL(Pxadj),  %VAL(Padjncy), %VAL(Pperm),
     & %VAL(Pinvp),  %VAL(Pcolcnt), nnzl, nsub, nsuper,
     & %VAL(Psnode), %VAL(Pxsuper), iwsiz, %VAL(Piwork), flag)
C
C CHECK ERROR FLAG
C
      if (flag .eq. -1) then
        call MEXERRMSGTXT('Insufficient working storage 
     &                     in sfinit.')
      endif
C
      Plindx  = MXCALLOC(2*nsub, 4)
      call symfct(
     & neqns,  anzmax,  %VAL(Pxadj),  %VAL(Padjncy), %VAL(Pperm), 
     & %VAL(Pinvp),%VAL(Pcolcnt),nsuper,%VAL(Pxsuper),%VAL(Psnode), 
     & nsub,  %VAL(Pxlindx),  %VAL(Plindx),  %VAL(Pxlnz), iwsiz, 
     & %VAL(Piwork),
     & flag )
C
C CHECK ERROR FLAG
C
      if (flag .eq. -1) then
        CALL MEXERRMSGTXT('Insufficient integer working space 
     &                     in symfct.')
      endif
      if (flag .eq. -2) then
        call MEXERRMSGTXT('Inconsistancy in the input in symfct.')
      endif
C
      cachsz = MXGETSCALAR(PRHS(4))
      call bfinit(neqns, nsuper, %VAL(Pxsuper),%VAL(Psnode),
     & %VAL(Pxlindx), %VAL(Plindx), cachsz, tmpsiz, %VAL(Psplit))
C
C CREATE MATRICES FOR RETURN ARGUMENTS
C
      PLHS(1)  = MXCREATEDOUBLEMATRIX(neqns+1, 1,0)
      PLHS(2)  = MXCREATEDOUBLEMATRIX(1,       1,0)
      PLHS(3)  = MXCREATEDOUBLEMATRIX(nsuper+1,1,0)
      PLHS(4)  = MXCREATEDOUBLEMATRIX(nsuper+1,1,0)
      PLHS(5)  = MXCREATEDOUBLEMATRIX(nsub,    1,0)
      PLHS(6)  = MXCREATEDOUBLEMATRIX(neqns,   1,0)
      PLHS(7)  = MXCREATEDOUBLEMATRIX(neqns,   1,0)
      PLHS(8)  = MXCREATEDOUBLEMATRIX(1,       1,0)
      PLHS(9)  = PRHS(2)
      PLHS(10) = PRHS(3)
C
C DEREFERENCE OUTPUT ARGUMENTS TO GET REAL PART POINTERS
C
      PMxlnz   = MXGETPR(PLHS(1))
      PMnnzl   = MXGETPR(PLHS(2))
      PMxsuper = MXGETPR(PLHS(3))
      PMxlindx = MXGETPR(PLHS(4))
      PMlindx  = MXGETPR(PLHS(5))
      PMsnode  = MXGETPR(PLHS(6))
      PMsplit  = MXGETPR(PLHS(7))
      PMtmpsiz = MXGETPR(PLHS(8))

      call int2real(neqns, nsuper, nsub,
     & %VAL( Pxlnz), nnzl, %VAL( Pxsuper),
     & %VAL( Pxlindx),%VAL( Plindx),%VAL( Psnode),%VAL( Psplit),
     & %VAL( Pperm),%VAL( Pinvp), tmpsiz,
     & %VAL(PMxlnz),%VAL(PMnnzl),%VAL(PMxsuper),
     & %VAL(PMxlindx),%VAL(PMlindx),%VAL(PMsnode),%VAL(PMsplit),
     & %VAL(PMperm),%VAL(PMinvp), %VAL(PMtmpsiz))

      RETURN
      END


c------------------------------------------------------
c Convert mxadj and madjncy of range [0:n-1]
c        to  xadj and  adjncy of range [1:n];
c Convert mperm and minvp of type real*8 in Matlab
c        to  perm and  invp of type integer in Fortran
c------------------------------------------------------
      subroutine convertin(neqns, anzmax,
     &           mxadj, madjncy, mperm, minvp,
     &            xadj,  adjncy,  perm,  invp)
      integer   neqns, anzmax, i
      integer   mxadj(neqns+1), madjncy(anzmax),
     &           xadj(neqns+1),  adjncy(anzmax)
      real*8    mperm(neqns), minvp(neqns)
      integer    perm(neqns),  invp(neqns)

      do 10 i = 1, neqns
         xadj(i) = mxadj(i) + 1
         perm(i) = mperm(i)
         invp(i) = minvp(i)
 10   continue
      xadj(neqns+1) = mxadj(neqns+1) + 1

      do 20 i = 1, anzmax
         adjncy(i) = madjncy(i) + 1
 20   continue
      return
      end

c-------------------------------------------------
c Convert outputs from integer type in Fortran
c                   to real*8  type in Matlab
c-------------------------------------------------
      subroutine int2real(neqns, nsuper, nsub,
     &            xlnz,  nnzl, xsuper, 
     &            xlindx, lindx, snode, split,
     &            perm,  invp, tmpsiz,
     &           mxlnz, mnnzl, mxsuper, 
     &           mxlindx,mlindx,msnode,msplit,
     &           mperm, minvp, mtmpsiz)

      integer  neqns, nsuper, nsub, i
      integer   xlnz(1),  nnzl,  tmpsiz,  xsuper(1) 
      real*8   mxlnz(1), mnnzl, mtmpsiz, mxsuper(1)
      integer   xlindx(1),  lindx(1),  snode(1),  split(1)
      real*8   mxlindx(1), mlindx(1), msnode(1), msplit(1)
      integer   perm(1),  invp(1)
      real*8   mperm(1), minvp(1)

      mnnzl = nnzl
      mtmpsiz = tmpsiz
      do 10 i = 1, neqns
         mxlnz(i)  = xlnz(i)
         msnode(i) = snode(i)
         msplit(i) = split(i)
         mperm(i)  = perm(i)
         minvp(i)  = invp(i)
 10   continue
      mxlnz(neqns+1) = xlnz(neqns+1)

      do 20 i = 1, nsuper+1
         mxsuper(i) = xsuper(i)
         mxlindx(i) = xlindx(i)
 20   continue

      do 30 i = 1, nsub
         mlindx(i) = lindx(i)
 30   continue
      return
      end
