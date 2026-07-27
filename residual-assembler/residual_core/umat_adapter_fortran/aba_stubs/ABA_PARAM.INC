C=======================================================================
C  ABA_PARAM.INC  -- standalone shim
C
C  Abaqus ships this include with its Fortran environment; it is NOT part
C  of the (MIT) UMAT sources, so we supply a compatible copy here so the
C  UMAT can be compiled outside Abaqus.  The two lines below are the
C  documented contents of the Abaqus double-precision ABA_PARAM.INC:
C    * IMPLICIT REAL*8(A-H,O-Z)  -> real*8 default for a-h, o-z names
C    * NPRECD = 2                -> "double precision" precision flag
C
C  Put this file's directory FIRST on the compiler include path (-I).
C=======================================================================
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (NPRECD=2)
