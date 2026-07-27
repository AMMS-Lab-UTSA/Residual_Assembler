C ======================================================================
C  elastic_oti_umat.for  (Milestone M4)
C
C  OTI-parameter-seeded isotropic linear-elastic verification UMAT. It is the
C  transformed counterpart of elastic_export_umat.for: instead of writing the
C  parameter derivatives dsigma/dE and dsigma/dnu ANALYTICALLY, it obtains them by
C  OTI automatic differentiation -- seeding E and nu as OTI numbers, running the
C  ordinary constitutive equations in OTI arithmetic, and extracting the first-
C  order coefficients. This is the mechanism that generalizes to crystal
C  plasticity (where analytic parameter derivatives are infeasible).
C
C  OTI library: OTIM4N1 (module OTIM4N1, type ONUMM4N1) = m=4 INDEPENDENT
C  first-order directions E1..E4. Seeding E along E1 and nu along E2 gives, in a
C  SINGLE constitutive evaluation:
C        SIG_OTI(I)%E1 = d sigma_i / dE
C        SIG_OTI(I)%E2 = d sigma_i / dnu
C  (%E3,%E4 are free for more parameters, e.g. C11/C12/C44 later.)
C
C  SDV contract (unchanged from M3):
C     STATEV(1:36)  = real DDSDDE, 6x6 row-major (= D_OTI%R, the exact real tangent)
C     STATEV(37:42) = dsigma/dE   (Voigt)
C     STATEV(43:48) = dsigma/dnu  (Voigt)
C
C  For this linear-elastic verification the OTI stress is formed from the TOTAL
C  end-of-increment strain eps = STRAN+DSTRAN (NOT an incremental update from the
C  incoming real STRESS, which carries no OTI coefficients from the previous
C  increment). The REAL stress is read back from the SAME OTI result so the real
C  and differentiated paths use identical equations.
C
C  PROPS(1) = E, PROPS(2) = nu.  *Depvar must be 48. Link with OTIM4N1 (see the
C  local abaqus_v6.env adding -I/-L/-lotim4n1). USE must precede ABA_PARAM.INC.
C ======================================================================
      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
C
      USE OTIM4N1
      INCLUDE 'ABA_PARAM.INC'
C
      CHARACTER*80 CMNAME
      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)
C
      PARAMETER (ZERO=0.D0, ONE=1.D0, TWO=2.D0)
      DIMENSION EPS(6)
      TYPE(ONUMM4N1) :: E_OTI, NU_OTI, LAM_OTI, MU_OTI
      TYPE(ONUMM4N1) :: D_OTI(6,6), SIG_OTI(6)
C
C  ---- seed E and nu SIMULTANEOUSLY along independent OTI directions --------
      E_OTI  = PROPS(1) + E1
      NU_OTI = PROPS(2) + E2
C
C  ---- OTI isotropic constants (engineering-shear convention) ---------------
      LAM_OTI = E_OTI*NU_OTI /
     &          ((ONE + NU_OTI)*(ONE - TWO*NU_OTI))
      MU_OTI  = E_OTI / (TWO*(ONE + NU_OTI))
C
C  ---- OTI tangent D_OTI (6x6) ----------------------------------------------
      DO I=1,NTENS
        DO J=1,NTENS
          D_OTI(I,J) = ZERO
        END DO
      END DO
      DO I=1,3
        DO J=1,3
          D_OTI(I,J) = LAM_OTI
        END DO
        D_OTI(I,I) = LAM_OTI + TWO*MU_OTI
      END DO
      DO I=4,NTENS
        D_OTI(I,I) = MU_OTI
      END DO
C
C  ---- total end-of-increment strain (real) ---------------------------------
      DO I=1,NTENS
        EPS(I) = STRAN(I) + DSTRAN(I)
      END DO
C
C  ---- OTI stress: sigma_OTI = D_OTI . eps -----------------------------------
      DO I=1,NTENS
        SIG_OTI(I) = ZERO
        DO J=1,NTENS
          SIG_OTI(I) = SIG_OTI(I) + D_OTI(I,J)*EPS(J)
        END DO
      END DO
C
C  ---- REAL stress from the SAME OTI result ---------------------------------
      DO I=1,NTENS
        STRESS(I) = SIG_OTI(I)%R
      END DO
C
C  ---- REAL DDSDDE (exact real tangent; unchanged from M3) -------------------
      DO I=1,NTENS
        DO J=1,NTENS
          DDSDDE(I,J) = D_OTI(I,J)%R
        END DO
      END DO
C
C  ---- pack SDVs ------------------------------------------------------------
C  STATEV(1:36) = DDSDDE row-major (index = 6*(i-1)+j)
      K = 0
      DO I=1,NTENS
        DO J=1,NTENS
          K = K + 1
          STATEV(K) = DDSDDE(I,J)
        END DO
      END DO
C  STATEV(37:42) = d sigma/dE (OTI direction 1), STATEV(43:48) = d sigma/dnu (dir 2)
      DO I=1,NTENS
        STATEV(36+I) = SIG_OTI(I)%E1
        STATEV(42+I) = SIG_OTI(I)%E2
      END DO
C
      RETURN
      END
