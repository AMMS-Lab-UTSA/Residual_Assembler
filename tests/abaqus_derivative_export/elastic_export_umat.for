C ======================================================================
C  elastic_export_umat.for
C
C  Small-strain isotropic linear-elastic UMAT used ONLY to verify the M3
C  ODB -> derivative-fields export pipeline (it is NOT the OTIS product).
C
C  It returns the exact stress and consistent tangent, and writes into the
C  state variables, at every integration point, the quantities the field-
C  driven residual method consumes:
C
C     STATEV(1:36)  = DDSDDE, 6x6 row-major   (i outer, j inner)
C     STATEV(37:42) = d sigma / dE            (Voigt, = sigma / E)
C     STATEV(43:48) = d sigma / d nu          (Voigt, = (dD/dnu) . eps)
C
C  Conventions match residual_core: Abaqus Voigt order (11,22,33,12,13,23),
C  ENGINEERING shear in STRAN/DSTRAN, so DDSDDE == core voigt.isotropic_D.
C
C  PROPS(1) = E,  PROPS(2) = nu.   *Depvar must be 48.
C ======================================================================
      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
C
      INCLUDE 'ABA_PARAM.INC'
C
      CHARACTER*80 CMNAME
      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)
C
      PARAMETER (ZERO=0.D0, ONE=1.D0, TWO=2.D0, THREE=3.D0)
      DIMENSION DDNU(6,6), EPS(6)
C
      EMOD = PROPS(1)
      ENU  = PROPS(2)
C
C  ---- isotropic elastic constants (engineering-shear convention) -----
      ELAM = EMOD*ENU/((ONE+ENU)*(ONE-TWO*ENU))
      EMU  = EMOD/(TWO*(ONE+ENU))
C
C  ---- consistent tangent DDSDDE = isotropic_D(E,nu) ------------------
      DO I=1,NTENS
        DO J=1,NTENS
          DDSDDE(I,J)=ZERO
        END DO
      END DO
      DO I=1,3
        DO J=1,3
          DDSDDE(I,J)=ELAM
        END DO
        DDSDDE(I,I)=ELAM+TWO*EMU
      END DO
      DO I=4,NTENS
        DDSDDE(I,I)=EMU
      END DO
C
C  ---- incremental stress update (linear elastic) --------------------
      DO I=1,NTENS
        DO J=1,NTENS
          STRESS(I)=STRESS(I)+DDSDDE(I,J)*DSTRAN(J)
        END DO
      END DO
C
C  ---- total strain at end of increment ------------------------------
      DO I=1,NTENS
        EPS(I)=STRAN(I)+DSTRAN(I)
      END DO
C
C  ---- d D / d nu  (analytic) ----------------------------------------
C  lam = E nu / ((1+nu)(1-2nu)),  mu = E / (2(1+nu))
C  d lam/d nu = E (1 + 2 nu^2) / ((1+nu)(1-2nu))^2
C  d mu /d nu = -E / (2 (1+nu)^2)
      DEN  = (ONE+ENU)*(ONE-TWO*ENU)
      DLAM = EMOD*(ONE+TWO*ENU*ENU)/(DEN*DEN)
      DMU  = -EMOD/(TWO*(ONE+ENU)*(ONE+ENU))
      DO I=1,NTENS
        DO J=1,NTENS
          DDNU(I,J)=ZERO
        END DO
      END DO
      DO I=1,3
        DO J=1,3
          DDNU(I,J)=DLAM
        END DO
        DDNU(I,I)=DLAM+TWO*DMU
      END DO
      DO I=4,NTENS
        DDNU(I,I)=DMU
      END DO
C
C  ---- write derivative fields into STATEV ---------------------------
C  STATEV(1:36) = DDSDDE row-major (index = 6*(i-1)+j)
      K=0
      DO I=1,NTENS
        DO J=1,NTENS
          K=K+1
          STATEV(K)=DDSDDE(I,J)
        END DO
      END DO
C  STATEV(37:42) = d sigma/dE = sigma / E   (D is linear in E)
      DO I=1,NTENS
        STATEV(36+I)=STRESS(I)/EMOD
      END DO
C  STATEV(43:48) = d sigma/d nu = (dD/dnu) . eps
      DO I=1,NTENS
        DSDNU=ZERO
        DO J=1,NTENS
          DSDNU=DSDNU+DDNU(I,J)*EPS(J)
        END DO
        STATEV(42+I)=DSDNU
      END DO
C
      RETURN
      END
