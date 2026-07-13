!=======================================================================
!  umat_mock.f90  --  a trivial *mock* UMAT with the EXACT Abaqus UMAT
!  argument list, used to exercise umat_driver + umat_replay.py end-to-end
!  WITHOUT ifort/Abaqus and WITHOUT the real crystal-plasticity source.
!
!  It is NOT the crystal-plasticity model.  It is small-strain isotropic
!  linear elasticity, wired to the SAME STATEV slots the real UMAT uses so
!  the replay/compare code path is identical:
!     STATEV(1:9)   = rotation matrix (identity here)
!     STATEV(35)    = a monotone "accumulated strain" scalar (history proof)
!     STATEV(48:53) = Cauchy stress  (matches umat.for: stress(K)=STATEV(47+K))
!     STATEV(81),(85),(89) = plastic F diagonal (identity here)
!
!  Small-strain measure:  eps = sym(F) - I,  stress = C_iso : eps
!  (engineering shear in Voigt, Abaqus order 11,22,33,12,13,23).
!  Identity deformation  ->  zero stress  (used by umat_replay --dry-run
!  as a known-answer sanity check).
!
!  Residual_Assembler, MIT.  Original work.
!=======================================================================
      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD, &
       RPL,DDSDDT,DRPLDE,DRPLDT, &
       STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME, &
       NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT, &
       CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
      implicit none
      character(len=80) :: CMNAME
      integer :: NDI,NSHR,NTENS,NSTATV,NPROPS,NOEL,NPT,LAYER,KSPT,KSTEP,KINC
      real(8) :: STRESS(NTENS),STATEV(NSTATV),DDSDDE(NTENS,NTENS)
      real(8) :: DDSDDT(NTENS),DRPLDE(NTENS),STRAN(NTENS),DSTRAN(NTENS)
      real(8) :: TIME(2),PREDEF(1),DPRED(1),PROPS(NPROPS)
      real(8) :: COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)
      real(8) :: SSE,SPD,SCD,RPL,DRPLDT,DTIME,TEMP,DTEMP,PNEWDT,CELENT

      real(8), parameter :: E = 1.0d5, nu = 0.3d0
      real(8) :: lam, mu, eps(6), deq
      integer :: i, j

      mu  = E / (2.0d0*(1.0d0+nu))
      lam = E*nu / ((1.0d0+nu)*(1.0d0-2.0d0*nu))

      ! ---- once-only initialisation (mirror umat.for kinc<=1 block) ----
      if (KINC <= 1 .and. KSTEP == 1) then
         do i = 1, NSTATV
            STATEV(i) = 0.0d0
         end do
         ! rotation matrix STATEV(1:9) from PROPS(2:10) if available
         if (NPROPS >= 10) then
            do i = 1, 3
               do j = 1, 3
                  STATEV(j+(i-1)*3) = PROPS(j+1+((i-1)*3))
               end do
            end do
         else
            STATEV(1) = 1.0d0; STATEV(5) = 1.0d0; STATEV(9) = 1.0d0
         end if
         STATEV(81) = 1.0d0; STATEV(85) = 1.0d0; STATEV(89) = 1.0d0
         STATEV(35) = 0.0d0
      end if

      ! ---- small strain from the total deformation gradient ----
      eps(1) = DFGRD1(1,1) - 1.0d0
      eps(2) = DFGRD1(2,2) - 1.0d0
      eps(3) = DFGRD1(3,3) - 1.0d0
      eps(4) = 0.5d0*(DFGRD1(1,2)+DFGRD1(2,1))   ! e12
      eps(5) = 0.5d0*(DFGRD1(1,3)+DFGRD1(3,1))   ! e13
      eps(6) = 0.5d0*(DFGRD1(2,3)+DFGRD1(3,2))   ! e23

      ! ---- isotropic tangent (engineering shear Voigt) ----
      do i = 1, NTENS
         do j = 1, NTENS
            DDSDDE(i,j) = 0.0d0
         end do
      end do
      do i = 1, 3
         do j = 1, 3
            DDSDDE(i,j) = lam
         end do
         DDSDDE(i,i) = lam + 2.0d0*mu
      end do
      if (NTENS >= 4) DDSDDE(4,4) = mu
      if (NTENS >= 5) DDSDDE(5,5) = mu
      if (NTENS >= 6) DDSDDE(6,6) = mu

      ! ---- stress = C : eps  (engineering shear -> use 2*e for shear) ----
      STRESS(1) = lam*(eps(1)+eps(2)+eps(3)) + 2.0d0*mu*eps(1)
      STRESS(2) = lam*(eps(1)+eps(2)+eps(3)) + 2.0d0*mu*eps(2)
      STRESS(3) = lam*(eps(1)+eps(2)+eps(3)) + 2.0d0*mu*eps(3)
      if (NTENS >= 4) STRESS(4) = 2.0d0*mu*eps(4)
      if (NTENS >= 5) STRESS(5) = 2.0d0*mu*eps(5)
      if (NTENS >= 6) STRESS(6) = 2.0d0*mu*eps(6)

      ! ---- write UMAT-compatible STATEV slots ----
      do i = 1, min(6,NTENS)
         STATEV(47+i) = STRESS(i)          ! Cauchy stress, real-UMAT slot
      end do
      deq = sqrt(eps(1)**2+eps(2)**2+eps(3)**2 &
               + 2.0d0*(eps(4)**2+eps(5)**2+eps(6)**2))
      STATEV(35) = STATEV(35) + deq*DTIME   ! monotone history witness

      PNEWDT = 1.0d0
      return
      END SUBROUTINE UMAT
