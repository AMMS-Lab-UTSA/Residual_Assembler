C ======================================================================
C  cubic_oriented_oti_umat.for  (Milestone M4 -- oriented FCC cubic crystal)
C
C  OTI-parameter-seeded CUBIC elasticity with a NON-TRIVIAL crystal ORIENTATION
C  (Bunge Euler angles phi1,PHI,phi2 = 25,40,15 deg). This closes the last cheap
C  gap before crystal plasticity: it proves crystal-orientation handling, not just
C  aligned [100] cubic arithmetic.
C
C  A SINGLE OTI-rotated operator is used for STRESS, DDSDDE and every dsigma/dCij:
C  the helper CUBIC_STRESS_OTI rotates the strain into the crystal frame
C  (eps_c = G eps G^T), applies the cubic law with the OTI-seeded constants, and
C  rotates the stress back (sig = G^T sig_c G). The real tangent DDSDDE is the
C  same operator's real part, built column-by-column from unit strains -- there is
C  NO separate real vs OTI rotation path.
C
C  SDV contract:
C     STATEV(1:36)  = real DDSDDE (rotated), 6x6 row-major
C     STATEV(37:42) = d sigma / dC11
C     STATEV(43:48) = d sigma / dC12
C     STATEV(49:54) = d sigma / dC44
C  PROPS(1)=C11, PROPS(2)=C12, PROPS(3)=C44.  *Depvar=54.  Link OTIM4N1.
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
      PARAMETER (PI=3.141592653589793D0)
      DIMENSION EPS(6), EV(6), G(3,3)
      TYPE(ONUMM4N1) :: C11O, C12O, C44O
      TYPE(ONUMM4N1) :: DG(6,6), SIG_OTI(6), COL(6)
C
C  ---- seed C11,C12,C44 along OTI directions E1,E2,E3 ------------------------
      C11O = PROPS(1) + E1
      C12O = PROPS(2) + E2
      C44O = PROPS(3) + E3
C
C  ---- fixed crystal orientation (Bunge phi1,PHI,phi2 = 25,40,15 deg) --------
      CALL MAKE_G(25.D0*PI/180.D0, 40.D0*PI/180.D0, 15.D0*PI/180.D0, G)
C
C  ---- build the rotated OTI stiffness DG(6,6) column-by-column --------------
C       column J = stress response to unit Voigt strain e_J (same operator).
      DO J=1,NTENS
        DO I=1,NTENS
          EV(I)=0.D0
        END DO
        EV(J)=1.D0
        CALL CUBIC_STRESS_OTI(EV, C11O, C12O, C44O, G, COL)
        DO I=1,NTENS
          DG(I,J)=COL(I)
        END DO
      END DO
C
C  ---- total end-of-increment strain (real, engineering shear) --------------
      DO I=1,NTENS
        EPS(I)=STRAN(I)+DSTRAN(I)
      END DO
C
C  ---- stress via the SAME operator:  sigma_OTI = DG . eps -------------------
      DO I=1,NTENS
        SIG_OTI(I)=0.D0
        DO J=1,NTENS
          SIG_OTI(I)=SIG_OTI(I)+DG(I,J)*EPS(J)
        END DO
      END DO
C
C  ---- real stress + real DDSDDE from the same operator ---------------------
      DO I=1,NTENS
        STRESS(I)=SIG_OTI(I)%R
        DO J=1,NTENS
          DDSDDE(I,J)=DG(I,J)%R
        END DO
      END DO
C
C  ---- pack SDVs ------------------------------------------------------------
      K=0
      DO I=1,NTENS
        DO J=1,NTENS
          K=K+1
          STATEV(K)=DDSDDE(I,J)
        END DO
      END DO
      DO I=1,NTENS
        STATEV(36+I)=SIG_OTI(I)%E1     ! d sigma / dC11
        STATEV(42+I)=SIG_OTI(I)%E2     ! d sigma / dC12
        STATEV(48+I)=SIG_OTI(I)%E3     ! d sigma / dC44
      END DO
C
      RETURN
      END
C ======================================================================
C  Bunge rotation g = Rz(phi2) Rx(PHI) Rz(phi1)  (sample -> crystal).
C  eps_crystal = g eps_sample g^T ,  sigma_sample = g^T sigma_crystal g.
C ======================================================================
      SUBROUTINE MAKE_G(P1, P, P2, G)
      INCLUDE 'ABA_PARAM.INC'
      DIMENSION G(3,3), RZ1(3,3), RXM(3,3), RZ2(3,3), T(3,3)
      CALL RZ(P1, RZ1)
      CALL RX(P,  RXM)
      CALL RZ(P2, RZ2)
      DO I=1,3
        DO J=1,3
          T(I,J)=0.D0
          DO KK=1,3
            T(I,J)=T(I,J)+RXM(I,KK)*RZ1(KK,J)
          END DO
        END DO
      END DO
      DO I=1,3
        DO J=1,3
          G(I,J)=0.D0
          DO KK=1,3
            G(I,J)=G(I,J)+RZ2(I,KK)*T(KK,J)
          END DO
        END DO
      END DO
      RETURN
      END

      SUBROUTINE RZ(A, R)
      INCLUDE 'ABA_PARAM.INC'
      DIMENSION R(3,3)
      C=COS(A)
      S=SIN(A)
      R(1,1)= C; R(1,2)= S; R(1,3)=0.D0
      R(2,1)=-S; R(2,2)= C; R(2,3)=0.D0
      R(3,1)=0.D0; R(3,2)=0.D0; R(3,3)=1.D0
      RETURN
      END

      SUBROUTINE RX(A, R)
      INCLUDE 'ABA_PARAM.INC'
      DIMENSION R(3,3)
      C=COS(A)
      S=SIN(A)
      R(1,1)=1.D0; R(1,2)=0.D0; R(1,3)=0.D0
      R(2,1)=0.D0; R(2,2)= C; R(2,3)= S
      R(3,1)=0.D0; R(3,2)=-S; R(3,3)= C
      RETURN
      END
C ======================================================================
C  CUBIC_STRESS_OTI: sigma_voigt(OTI) from strain_voigt(real, engineering shear)
C  eps_c = G eps G^T (real) -> cubic law (OTI) -> sigma = G^T sig_c G (OTI).
C ======================================================================
      SUBROUTINE CUBIC_STRESS_OTI(EV, C11O, C12O, C44O, G, SV)
      USE OTIM4N1
      INCLUDE 'ABA_PARAM.INC'
      DIMENSION EV(6), G(3,3), ET(3,3), EC(3,3), TMP(3,3)
      TYPE(ONUMM4N1) :: C11O, C12O, C44O
      TYPE(ONUMM4N1) :: SV(6), SC(3,3), SS(3,3), TO(3,3)
C  strain Voigt (engineering shear) -> tensor
      ET(1,1)=EV(1); ET(2,2)=EV(2); ET(3,3)=EV(3)
      ET(1,2)=0.5D0*EV(4); ET(2,1)=ET(1,2)
      ET(1,3)=0.5D0*EV(5); ET(3,1)=ET(1,3)
      ET(2,3)=0.5D0*EV(6); ET(3,2)=ET(2,3)
C  eps_c = G eps G^T  (all real)
      DO I=1,3
        DO J=1,3
          TMP(I,J)=0.D0
          DO KK=1,3
            TMP(I,J)=TMP(I,J)+G(I,KK)*ET(KK,J)
          END DO
        END DO
      END DO
      DO I=1,3
        DO J=1,3
          EC(I,J)=0.D0
          DO KK=1,3
            EC(I,J)=EC(I,J)+TMP(I,KK)*G(J,KK)
          END DO
        END DO
      END DO
C  cubic law in the crystal frame (OTI)
      SC(1,1)=C11O*EC(1,1)+C12O*EC(2,2)+C12O*EC(3,3)
      SC(2,2)=C12O*EC(1,1)+C11O*EC(2,2)+C12O*EC(3,3)
      SC(3,3)=C12O*EC(1,1)+C12O*EC(2,2)+C11O*EC(3,3)
      SC(1,2)=(2.D0*C44O)*EC(1,2); SC(2,1)=SC(1,2)
      SC(1,3)=(2.D0*C44O)*EC(1,3); SC(3,1)=SC(1,3)
      SC(2,3)=(2.D0*C44O)*EC(2,3); SC(3,2)=SC(2,3)
C  sigma_sample = G^T sig_c G  (OTI; G real)
      DO I=1,3
        DO J=1,3
          TO(I,J)=0.D0
          DO KK=1,3
            TO(I,J)=TO(I,J)+G(KK,I)*SC(KK,J)
          END DO
        END DO
      END DO
      DO I=1,3
        DO J=1,3
          SS(I,J)=0.D0
          DO KK=1,3
            SS(I,J)=SS(I,J)+TO(I,KK)*G(KK,J)
          END DO
        END DO
      END DO
C  stress tensor -> Voigt
      SV(1)=SS(1,1); SV(2)=SS(2,2); SV(3)=SS(3,3)
      SV(4)=SS(1,2); SV(5)=SS(1,3); SV(6)=SS(2,3)
      RETURN
      END
