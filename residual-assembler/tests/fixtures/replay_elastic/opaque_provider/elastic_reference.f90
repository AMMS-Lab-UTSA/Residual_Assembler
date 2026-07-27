! =====================================================================
!  elastic_reference.f90
!  REFERENCE ABI provider for the resasm replay tool  --  NOT the real
!  JHU OTI binary.  It implements the resasm_mat_abi.h v1 C ABI for an
!  isotropic linear-elastic small-strain material, computing the parameter
!  derivatives d(stress)/d(E,nu) ANALYTICALLY (closed form) rather than by
!  OTI arithmetic.
!
!  Its ONLY purpose is to exercise the collaborator-side residual/sensitivity
!  pipeline (adapter -> replay -> R,K,R_p -> du/dp -> dq/dp) through the exact
!  compiled C boundary the real binary will use, so the real OTI object can
!  replace this file WITHOUT any change to the Python residual code.
!
!  It deliberately keeps NO OTI type visible: it exposes real stress and, in a
!  SEPARATE array, the first-order derivative coefficients -- exactly as the
!  ABI requires.  Voigt order [11,22,33,12,13,23], engineering shear, matching
!  residual_core/formulations/c3d8_kernel.isotropic_D.
! =====================================================================
module resasm_elastic_reference
  use, intrinsic :: iso_c_binding
  implicit none

  integer(c_int), parameter :: ABI_VERSION = 1
  integer(c_int), parameter :: KIN_SMALL   = 0
  integer(c_int), parameter :: OK = 0, ERR_ABI = 1, ERR_DIMS = 2, ERR_NULL = 3, &
                               ERR_SEED = 4, ERR_KIN = 5, ERR_STATE = 6

  type, bind(C) :: resasm_mat_desc_t
     integer(c_int) :: abi_version, ntens, nprops, nstatev, kinematics, order
  end type resasm_mat_desc_t

contains

  ! Lame parameters and their derivatives wrt E (which=1) or nu (which=2).
  subroutine lame_and_deriv(E, nu, which, lam, mu, dlam, dmu)
    real(c_double), intent(in)  :: E, nu
    integer,        intent(in)  :: which
    real(c_double), intent(out) :: lam, mu, dlam, dmu
    real(c_double) :: den
    den = (1.0d0 + nu) * (1.0d0 - 2.0d0 * nu)
    lam = E * nu / den
    mu  = E / (2.0d0 * (1.0d0 + nu))
    if (which == 1) then                 ! d/dE  (both linear in E)
       dlam = lam / E
       dmu  = mu / E
    else                                 ! d/dnu
       dlam = E * (1.0d0 + 2.0d0 * nu * nu) / (den * den)
       dmu  = -E / (2.0d0 * (1.0d0 + nu) * (1.0d0 + nu))
    end if
  end subroutine lame_and_deriv

  ! Assemble the 6x6 elastic tangent (or its parameter derivative) from lam,mu.
  subroutine build_D(lam, mu, D)
    real(c_double), intent(in)  :: lam, mu
    real(c_double), intent(out) :: D(6,6)
    integer :: i, j
    D = 0.0d0
    do i = 1, 3
       do j = 1, 3
          D(i,j) = lam
       end do
       D(i,i) = lam + 2.0d0 * mu
    end do
    do i = 4, 6
       D(i,i) = mu
    end do
  end subroutine build_D

  ! ---- mat_describe_v1 -------------------------------------------------
  function mat_describe_v1(desc_out, model_id, model_id_cap) result(rc) &
       bind(C, name="mat_describe_v1")
    type(resasm_mat_desc_t), intent(out) :: desc_out
    character(kind=c_char), intent(out)  :: model_id(*)
    integer(c_int), value                :: model_id_cap
    integer(c_int) :: rc
    character(len=*), parameter :: name = "reference_elastic_isotropic"
    integer :: i, n
    desc_out%abi_version = ABI_VERSION
    desc_out%ntens = 6
    desc_out%nprops = 2
    desc_out%nstatev = 0
    desc_out%kinematics = KIN_SMALL
    desc_out%order = 1
    if (model_id_cap > 0) then
       n = min(len(name), int(model_id_cap) - 1)
       do i = 1, n
          model_id(i) = name(i:i)
       end do
       model_id(n+1) = c_null_char
    end if
    rc = OK
  end function mat_describe_v1

  ! ---- mat_eval_v1 -----------------------------------------------------
  function mat_eval_v1(desc, props, seed_indices, nseed, kin, dkin_dseed, &
       state_in, dstate_dseed_in, time_data, stress, dstress_dseed, &
       state_out, dstate_dseed_out, ddsdde, status) result(rc) &
       bind(C, name="mat_eval_v1")
    type(resasm_mat_desc_t), intent(in) :: desc
    real(c_double), intent(in)  :: props(*)
    integer(c_int), intent(in)  :: seed_indices(*)
    integer(c_int), value       :: nseed
    real(c_double), intent(in)  :: kin(*)
    type(c_ptr),    value       :: dkin_dseed        ! optional (unused: elastic)
    type(c_ptr),    value       :: state_in          ! optional (nstatev=0)
    type(c_ptr),    value       :: dstate_dseed_in   ! optional
    real(c_double), intent(in)  :: time_data(*)
    real(c_double), intent(out) :: stress(*)
    real(c_double), intent(out) :: dstress_dseed(*)
    type(c_ptr),    value       :: state_out         ! optional
    type(c_ptr),    value       :: dstate_dseed_out  ! optional
    real(c_double), intent(out) :: ddsdde(*)
    type(c_ptr),    value       :: status            ! optional
    integer(c_int) :: rc

    real(c_double) :: E, nu, lam, mu, dlam, dmu
    real(c_double) :: D(6,6), dD(6,6), eps(6)
    integer :: i, j, s, pidx
    integer(c_int), pointer :: pstatus

    rc = OK
    ! --- validate dimensions / contract --------------------------------
    ! This reference implements EXACTLY two parameters (E=props(1), nu=props(2));
    ! reject nprops/=2 so a seed index >2 cannot be silently mapped to d/dnu.
    if (desc%abi_version /= ABI_VERSION) then; rc = ERR_ABI;  call finish(); return; end if
    if (desc%ntens /= 6 .or. desc%nprops /= 2) then; rc = ERR_DIMS; call finish(); return; end if
    if (desc%kinematics /= KIN_SMALL) then; rc = ERR_KIN; call finish(); return; end if
    if (desc%nstatev /= 0) then; rc = ERR_STATE; call finish(); return; end if
    if (nseed < 1) then; rc = ERR_SEED; call finish(); return; end if
    do s = 1, nseed
       if (seed_indices(s) < 1 .or. seed_indices(s) > desc%nprops) then
          rc = ERR_SEED; call finish(); return
       end if
    end do

    E  = props(1)
    nu = props(2)
    if (E <= 0.0d0 .or. nu <= -1.0d0 .or. nu >= 0.5d0) then
       rc = ERR_DIMS; call finish(); return
    end if

    ! --- real response: sigma = D(E,nu) . eps , DDSDDE = D --------------
    call lame_and_deriv(E, nu, 1, lam, mu, dlam, dmu)   ! lam,mu (derivs unused here)
    call build_D(lam, mu, D)
    eps = kin(1:6)
    do i = 1, 6
       stress(i) = 0.0d0
       do j = 1, 6
          stress(i) = stress(i) + D(i,j) * eps(j)
       end do
    end do
    do i = 1, 6
       do j = 1, 6
          ddsdde((i-1)*6 + j) = D(i,j)     ! row-major, per the ABI
       end do
    end do

    ! --- derivative coefficients: dsigma/dp_s = (dD/dp_s) . eps ---------
    do s = 1, nseed
       pidx = seed_indices(s)              ! 1 -> E, 2 -> nu (this reference)
       call lame_and_deriv(E, nu, pidx, lam, mu, dlam, dmu)
       call build_D(dlam, dmu, dD)         ! dD/dp has the same structure in (dlam,dmu)
       do i = 1, 6
          dstress_dseed((i-1)*nseed + s) = 0.0d0
          do j = 1, 6
             dstress_dseed((i-1)*nseed + s) = dstress_dseed((i-1)*nseed + s) + dD(i,j) * eps(j)
          end do
       end do
    end do

    ! elastic: no state to update; leave optional out-buffers untouched.
    call finish()
  contains
    subroutine finish()
      if (c_associated(status)) then
         call c_f_pointer(status, pstatus)
         pstatus = rc
      end if
    end subroutine finish
  end function mat_eval_v1

end module resasm_elastic_reference
