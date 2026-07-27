C ======================================================================
C  cubic_oti_umat.for  (Milestone M4 -- FCC / cubic single crystal)
C
C  OTI-parameter-seeded CUBIC (anisotropic) linear-elastic UMAT. Crystal axes are
C  aligned with the global axes ([100] simple tension), so the Voigt stiffness in
C  the crystal frame is used directly. This is the C11/C12/C44 configuration for
C  the FCC simpleTension case: it seeds the three independent cubic elastic
C  constants SIMULTANEOUSLY with OTIM4N1 directions E1, E2, E3 and extracts their
C  stress derivatives by automatic differentiation.
C
C  Cubic stiffness (Abaqus Voigt 11,22,33,12,13,23; engineering shear):
C        [ C11 C12 C12              ]
C        [ C12 C11 C12              ]
C        [ C12 C12 C11              ]
C        [             C44          ]
C        [                 C44      ]
C        [                     C44  ]
C
C  SDV contract:
C     STATEV(1:36)  = real DDSDDE, 6x6 row-major (= D_OTI%R)
C     STATEV(37:42) = d sigma / dC11   (Voigt)
C     STATEV(43:48) = d sigma / dC12   (Voigt)
C     STATEV(49:54) = d sigma / dC44   (Voigt)
C
C  PROPS(1)=C11, PROPS(2)=C12, PROPS(3)=C44.  *Depvar must be 54.
C  Link with OTIM4N1 (local abaqus_v6.env: -I/-L/-lotim4n1). USE precedes
C  ABA_PARAM.INC. Real stress read from the same OTI result (STRESS=SIG_OTI%R);
C  OTI stress formed from the total end-of-increment strain eps=STRAN+DSTRAN.
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
      PARAMETER (ZERO=0.D0)
      DIMENSION EPS(6)
      TYPE(ONUMM4N1) :: C11O, C12O, C44O
      TYPE(ONUMM4N1) :: D_OTI(6,6), SIG_OTI(6)
C
C  ---- seed C11, C12, C44 SIMULTANEOUSLY along independent OTI directions ----
      C11O = PROPS(1) + E1
      C12O = PROPS(2) + E2
      C44O = PROPS(3) + E3
C
C  ---- OTI cubic tangent D_OTI (6x6) ----------------------------------------
      DO I=1,NTENS
        DO J=1,NTENS
          D_OTI(I,J) = ZERO
        END DO
      END DO
      DO I=1,3
        DO J=1,3
          D_OTI(I,J) = C12O
        END DO
        D_OTI(I,I) = C11O
      END DO
      DO I=4,NTENS
        D_OTI(I,I) = C44O
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
C  ---- REAL stress + REAL DDSDDE from the SAME OTI result -------------------
      DO I=1,NTENS
        STRESS(I) = SIG_OTI(I)%R
        DO J=1,NTENS
          DDSDDE(I,J) = D_OTI(I,J)%R
        END DO
      END DO
C
C  ---- pack SDVs ------------------------------------------------------------
      K = 0
      DO I=1,NTENS
        DO J=1,NTENS
          K = K + 1
          STATEV(K) = DDSDDE(I,J)
        END DO
      END DO
      DO I=1,NTENS
        STATEV(36+I) = SIG_OTI(I)%E1      ! d sigma / dC11
        STATEV(42+I) = SIG_OTI(I)%E2      ! d sigma / dC12
        STATEV(48+I) = SIG_OTI(I)%E3      ! d sigma / dC44
      END DO
C
      RETURN
      END
