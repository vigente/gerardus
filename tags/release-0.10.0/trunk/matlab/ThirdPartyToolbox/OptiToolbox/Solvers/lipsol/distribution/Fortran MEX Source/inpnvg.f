C***********************************************************************
C
C   Version:        Beta-2.2
C   Last modified:  January 13, 1995
C   Author:         Yin Zhang
C                   Department of Mathematics and Statistics
C                   University of Maryland Baltimore County
C
C***********************************************************************
C Usage: [lnz] = inpnv(Pdiag,P,invp,perm,xlnz,xsuper,xlindx,lindx,nnzl)
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
C KEEP THE ABOVALE SUBROUTINE, ARGUMENT, AND FUNCTION DECLARATIONS 
C FOR USE IN ALL YOUR FORTRAN MEX FILES.
C---------------------------------------------------------------------
C
      INTEGER*4 PPdiag, PP, PMip, PMjp, Pip, Pjp
      INTEGER*4 PMinvp, PMperm,PMxlnz, PMxsuper
      INTEGER*4  Pinvp,  Pperm, Pxlnz,  Pxsuper
      INTEGER*4 PMxlindx,PMlindx
      INTEGER*4  Pxlindx, Plindx
      INTEGER*4 Plnz, Plink
      INTEGER*4 Pxadjf, Padjf, Panzf
      INTEGER neqns, nnzl, nsuper,maxsub, nzmax
C
C CHECK FOR PROPER NUMBER OF ARGUMENTS
C
      IF (NRHS .NE. 9) THEN
        CALL MEXERRMSGTXT('INPNV requires 9 input arguments')
      ELSEIF (NLHS .NE. 1) THEN
        CALL MEXERRMSGTXT('INPNV requires 1 output arguments')
      ENDIF
C
C SPECIFY THE DIMENSIION OF WORKING VECTORS
C
      neqns  = max(MXGETM(PRHS(2)),MXGETN(PRHS(2)))
      nsuper = max(MXGETM(PRHS(6)),MXGETN(PRHS(6))) - 1
      maxsub = max(MXGETM(PRHS(8)),MXGETN(PRHS(8)))
      nnzl   = MXGETSCALAR(PRHS(9))
      nzmax  = MXGETNZMAX(PRHS(2))
C
C CREATE OUTPUT PARAMETERS
C
      PLHS(1) = MXCREATEDOUBLEMATRIX(nnzl,1,0)
C
C CREATE WORKING PARAMETERS
C
      Pip     = MXCALLOC(nzmax,    4)
      Pjp     = MXCALLOC(neqns+1,  4)
      Pinvp   = MXCALLOC(neqns,    4)
      Pperm   = MXCALLOC(neqns,    4)
      Pxlnz   = MXCALLOC(neqns+1,  4)
      Pxsuper = MXCALLOC(nsuper+1, 4)
      Pxlindx = MXCALLOC(nsuper+1, 4)
      Plindx  = MXCALLOC(maxsub ,  4)
      Pxadjf  = MXCALLOC(neqns+1,  4)
      Plink   = MXCALLOC(neqns,    4)
      Padjf   = MXCALLOC(nzmax+neqns, 4)
      Panzf   = MXCALLOC(nzmax+neqns, 8)
C
C DEREFERENCE ARGUMENTS TO GET ARRAY POINTERS
C
      PPdiag   = MXGETPR(PRHS(1))
      PP       = MXGETPR(PRHS(2))
      PMip     = MXGETIR(PRHS(2))
      PMjp     = MXGETJC(PRHS(2))
      PMinvp   = MXGETPR(PRHS(3))
      PMperm   = MXGETPR(PRHS(4))
      PMxlnz   = MXGETPR(PRHS(5))
      PMxsuper = MXGETPR(PRHS(6))
      PMxlindx = MXGETPR(PRHS(7))
      PMlindx  = MXGETPR(PRHS(8))
      Plnz     = MXGETPR(PLHS(1))
C
C INPUT DATA TRANSFORMATION
C
       call convertin(neqns, maxsub, nzmax, nsuper,
     & %VAL(PMip),  %VAL(PMjp), %VAL(PMinvp), %VAL(PMperm),
     & %VAL( Pip),  %VAL( Pjp), %VAL( Pinvp), %VAL( Pperm),
     & %VAL(PMxlnz),%VAL(PMxsuper),%VAL(PMxlindx), %VAL(PMlindx), 
     & %VAL( Pxlnz),%VAL( Pxsuper),%VAL( Pxlindx), %VAL( Plindx))
C
C DO THE ACTUAL COMPUTATIONS IN A FORTRAN SUBROUTINE
C
      call addiag(neqns,
     &     %VAL(Pjp), %VAL(Pip), %VAL(PP), %VAL(PPdiag), 
     &     %VAL(Pxadjf), %VAL(Padjf), %VAL(Panzf))

      call inpnv(neqns,   %VAL(Pxadjf),  %VAL(Padjf), 
     &     %VAL(Panzf),   %VAL(Pperm),   %VAL(Pinvp), nsuper,
     &     %VAL(Pxsuper), %VAL(Pxlindx), %VAL(Plindx), 
     &     %VAL(Pxlnz),   %VAL(Plnz),    %VAL(Plink))

      RETURN
      END

c----------------------------------------------------------------
c Convert from C range [0:n-1] to Fortran range [1:n];
c Convert from type real*8 in Matlab to type integer in Fortran
c----------------------------------------------------------------
      subroutine convertin(neqns, maxsub, nzmax, nsuper,
     &           mip, mjp, minvp, mperm,  
     &            ip,  jp,  invp,  perm,  
     &           mxlnz, mxsuper, mxlindx,   mlindx,
     &            xlnz,  xsuper,  xlindx,    lindx)

      integer neqns, maxsub, nzmax, nsuper, i
      integer mip(nzmax), mjp(neqns+1), ip(nzmax), jp(neqns+1)

      real*8  minvp(neqns), mperm(neqns), mxlnz(neqns+1)
      integer  invp(neqns),  perm(neqns),  xlnz(neqns+1)

      real*8  mxsuper(nsuper+1),mxlindx(neqns+1),mlindx(maxsub)
      integer  xsuper(nsuper+1), xlindx(neqns+1), lindx(maxsub)

      do 10 i = 1, neqns
         jp(i)     = mjp(i) + 1
         invp(i)   = minvp(i)
         perm(i)   = mperm(i)
         xlnz(i)   = mxlnz(i)
 10   continue
      jp(neqns+1)     = mjp(neqns+1) + 1
      xlnz(neqns+1)   = mxlnz(neqns+1)

      do 15 i = 1, nsuper+1
         xsuper(i) = mxsuper(i)
         xlindx(i)   = mxlindx(i)
 15   continue

      do 20 i = 1, nzmax
         ip(i) = mip(i) + 1
 20   continue

      do 30 i = 1, maxsub
         lindx(i) = mlindx(i)
 30   continue
      return
      end
c***********************************************************************
c
c   Last modified:  November, 1994
c   Author:         Detong Zhang
c                   Department of Mathematics and Statistics
c                   University of Maryland Baltimore County
c
c***********************************************************************
      subroutine addiag(n,jp,ip,P,Pdiag,xadjf,adjf,anzf)
      integer n,jp(n+1),xadjf(n+1),ip(*),adjf(*)
      real*8  P(*), Pdiag(n), anzf(*)
c
      do 40 j=1,n
         anzf(jp(j)+j-1)=Pdiag(j)
         adjf(jp(j)+j-1)=j
         do 50 ij=jp(j),jp(j+1)-1
            anzf(ij+j)=P(ij)
            adjf(ij+j)=ip(ij)
   50    continue
         xadjf(j)=jp(j)+j-1
   40  continue
       xadjf(n+1)=jp(n+1)+n
       end
