! Driver for claim 2 (slides 26-27): OTI vs hand-coded analytical Jacobian of
! the crystal-plasticity flow rule computeFlowRule (12 slip systems), and
! centred / forward finite differences of the ANALYTICAL routine's primal output.
!
! It is compiled together with the flow-rule sources (not part of this
! repository; see presentation/claim2_flowrule_jacobian.py) with ONE compiler
! and ONE flag set for every method.
!
!   ./claim2 values      -> CSV of analytical, OTI and FD derivatives (stdout)
!   ./claim2 timing N R  -> per-call seconds of each method, R repetitions of N calls
program claim2_flowrule
  use Kinds
  use FlowRule_OTIS
  implicit none
  integer, parameter :: ns = 12, nv = 5, nsteps = 7
  real(kind=wp), parameter :: steps(nsteps) = (/ 1.0e-2_wp, 1.0e-3_wp, 1.0e-4_wp, 1.0e-5_wp, &
                                                 1.0e-6_wp, 1.0e-7_wp, 1.0e-8_wp /)
  real(kind=wp) :: tau(ns), temp, cdeg, back(ns), tssd(ns), tcs(ns), adi(ns), assd(ns), tgnd(ns), agnd(ns), dtime
  real(kind=wp) :: inc(ns), rtdi(ns), d1(ns), d2(ns), d3(ns), d4(ns), d5(ns), o1(ns), o2(ns), o3(ns)
  real(kind=wp) :: ana(ns, nv), oti(ns, nv), incA(ns), incO(ns)
  real(kind=wp) :: z0(ns, nv), z(ns, nv), fp(ns), fm(ns), f0(ns), h(ns), scale(ns), checksum
  real(kind=wp) :: t0, t1
  character(len=16) :: mode, arg
  integer :: v, k, a, iter, rep, ncalls, nrep
  integer(8) :: c0, c1, rate

  call get_command_argument(1, mode)
  tau = (/ 351341296.634225_wp,-18069556.4068816_wp,-333271740.227343_wp, &
           353247960.639986_wp,-16258249.2359092_wp,-336989711.404077_wp, &
           348287762.056055_wp, 3723893.08626250_wp,-352011655.142318_wp, &
           353955103.276707_wp, 1921201.35165965_wp,-355876304.628367_wp /)
  temp = 1123.0_wp; cdeg = 1.0_wp; back = 0.0_wp
  tssd = (/ 539422634.907665_wp, 539422634.907665_wp, 539422634.907665_wp, &
            512710648.819887_wp, 512710648.819887_wp, 512710648.819887_wp, &
            514299352.207645_wp, 514299352.207645_wp, 514299352.207645_wp, &
            538152337.407513_wp, 538152337.407513_wp, 538152337.407513_wp /)
  tcs = 297692600.395148_wp; adi = 0.0_wp; assd = 0.0_wp; tgnd = 0.0_wp; agnd = 0.0_wp
  dtime = 1.0_wp
  checksum = 0.0_wp

  ! variables differentiated, in the OTI seeding order E1..E5
  z0(:, 1) = tau; z0(:, 2) = back; z0(:, 3) = tssd; z0(:, 4) = tcs; z0(:, 5) = assd
  ! finite-difference scale: the variable itself, or the resolved shear stress
  ! of that slip system where the variable is zero
  scale = abs(tau)

  if (trim(mode) == 'values') then
    call analytical(z0, .true., incA, ana)
    call otis(z0, incO, oti)
    write(*, '(A)') 'kind,variable,step,slip_system,value'
    do a = 1, ns
      write(*, '(A,",",I1,",",ES10.3,",",I2,",",ES25.17)') 'primal_analytical', 0, 0.0_wp, a, incA(a)
      write(*, '(A,",",I1,",",ES10.3,",",I2,",",ES25.17)') 'primal_oti', 0, 0.0_wp, a, incO(a)
    end do
    do v = 1, nv
      do a = 1, ns
        write(*, '(A,",",I1,",",ES10.3,",",I2,",",ES25.17)') 'analytical', v, 0.0_wp, a, ana(a, v)
        write(*, '(A,",",I1,",",ES10.3,",",I2,",",ES25.17)') 'oti', v, 0.0_wp, a, oti(a, v)
      end do
      do k = 1, nsteps
        h = steps(k) * max(abs(z0(:, v)), scale)
        ! one slip system at a time: no assumption that the systems decouple
        do a = 1, ns
          z = z0; z(a, v) = z0(a, v) + h(a); call primal(z, fp)
          z = z0; z(a, v) = z0(a, v) - h(a); call primal(z, fm)
          call primal(z0, f0)
          write(*, '(A,",",I1,",",ES10.3,",",I2,",",ES25.17)') 'fd_centred', v, steps(k), a, &
               (fp(a) - fm(a)) / (2.0_wp * h(a))
          write(*, '(A,",",I1,",",ES10.3,",",I2,",",ES25.17)') 'fd_forward', v, steps(k), a, &
               (fp(a) - f0(a)) / h(a)
          write(*, '(A,",",I1,",",ES10.3,",",I2,",",ES25.17)') 'fd_backward', v, steps(k), a, &
               (f0(a) - fm(a)) / h(a)
        end do
        ! all slip systems at once (what a cheap FD Jacobian would do)
        z = z0; z(:, v) = z0(:, v) + h; call primal(z, fp)
        z = z0; z(:, v) = z0(:, v) - h; call primal(z, fm)
        do a = 1, ns
          write(*, '(A,",",I1,",",ES10.3,",",I2,",",ES25.17)') 'fd_centred_all_systems', v, &
               steps(k), a, (fp(a) - fm(a)) / (2.0_wp * h(a))
        end do
      end do
    end do

  else if (trim(mode) == 'timing') then
    call get_command_argument(2, arg); read(arg, *) ncalls
    call get_command_argument(3, arg); read(arg, *) nrep
    k = 4                                       ! the FD step used for "known step" timing
    do rep = 1, nrep
      ! analytical Jacobian (the hand-coded one), one call
      call analytical(z0, .true., incA, ana)                      ! warm-up
      call system_clock(c0, rate)
      do iter = 1, ncalls
        call analytical(z0, .true., incA, ana); checksum = checksum + ana(1, 1)
      end do
      call system_clock(c1)
      write(*, '(A,",",I3,",",ES24.16)') 'analytical', rep, dble(c1 - c0) / dble(rate) / dble(ncalls)
      ! OTI, one call
      call otis(z0, incO, oti)
      call system_clock(c0, rate)
      do iter = 1, ncalls
        call otis(z0, incO, oti); checksum = checksum + oti(1, 1)
      end do
      call system_clock(c1)
      write(*, '(A,",",I3,",",ES24.16)') 'oti', rep, dble(c1 - c0) / dble(rate) / dble(ncalls)
      ! primal only (one evaluation of the analytical routine without Jacobian)
      call system_clock(c0, rate)
      do iter = 1, ncalls
        call primal(z0, f0); checksum = checksum + f0(1)
      end do
      call system_clock(c1)
      write(*, '(A,",",I3,",",ES24.16)') 'primal', rep, dble(c1 - c0) / dble(rate) / dble(ncalls)
      ! centred FD Jacobian with a known step: 2 * nv primal calls (all systems at once)
      call system_clock(c0, rate)
      do iter = 1, ncalls
        do v = 1, nv
          h = steps(k) * max(abs(z0(:, v)), scale)
          z = z0; z(:, v) = z0(:, v) + h; call primal(z, fp)
          z = z0; z(:, v) = z0(:, v) - h; call primal(z, fm)
          checksum = checksum + (fp(1) - fm(1)) / (2.0_wp * h(1))
        end do
      end do
      call system_clock(c1)
      write(*, '(A,",",I3,",",ES24.16)') 'fd_centred_known_step', rep, dble(c1 - c0) / dble(rate) / dble(ncalls)
      ! forward FD Jacobian with a known step: nv + 1 primal calls
      call system_clock(c0, rate)
      do iter = 1, ncalls
        call primal(z0, f0)
        do v = 1, nv
          h = steps(k) * max(abs(z0(:, v)), scale)
          z = z0; z(:, v) = z0(:, v) + h; call primal(z, fp)
          checksum = checksum + (fp(1) - f0(1)) / h(1)
        end do
      end do
      call system_clock(c1)
      write(*, '(A,",",I3,",",ES24.16)') 'fd_forward_known_step', rep, dble(c1 - c0) / dble(rate) / dble(ncalls)
      ! centred FD including the step search over the whole ladder
      call system_clock(c0, rate)
      do iter = 1, ncalls / 10
        do k = 1, nsteps
          do v = 1, nv
            h = steps(k) * max(abs(z0(:, v)), scale)
            z = z0; z(:, v) = z0(:, v) + h; call primal(z, fp)
            z = z0; z(:, v) = z0(:, v) - h; call primal(z, fm)
            checksum = checksum + (fp(1) - fm(1)) / (2.0_wp * h(1))
          end do
        end do
      end do
      call system_clock(c1)
      k = 4
      write(*, '(A,",",I3,",",ES24.16)') 'fd_centred_step_search', rep, &
           dble(c1 - c0) / dble(rate) / dble(ncalls / 10)
    end do
    call cpu_time(t0); call cpu_time(t1)
    write(*, '(A,",",I3,",",ES24.16)') 'checksum', 0, checksum
  else
    write(*, '(A)') 'usage: claim2 values | claim2 timing NCALLS NREP'
    stop 2
  end if

contains

  subroutine unpack(z)
    real(kind=wp), intent(in) :: z(ns, nv)
    tau = z(:, 1); back = z(:, 2); tssd = z(:, 3); tcs = z(:, 4); assd = z(:, 5)
  end subroutine unpack

  subroutine analytical(z, jac, incr, der)
    real(kind=wp), intent(in) :: z(ns, nv)
    logical, intent(in) :: jac
    real(kind=wp), intent(out) :: incr(ns), der(ns, nv)
    call unpack(z)
    call computeFlowRule(tau, temp, cdeg, back, tssd, tcs, adi, assd, tgnd, agnd, dtime, &
                         incr, rtdi, d1, d2, d3, d4, d5, jac, o1, o2, o3)
    der(:, 1) = d1; der(:, 2) = d2; der(:, 3) = d3; der(:, 4) = d4; der(:, 5) = d5
  end subroutine analytical

  subroutine primal(z, incr)
    real(kind=wp), intent(in) :: z(ns, nv)
    real(kind=wp), intent(out) :: incr(ns)
    call unpack(z)
    call computeFlowRule(tau, temp, cdeg, back, tssd, tcs, adi, assd, tgnd, agnd, dtime, &
                         incr, rtdi, d1, d2, d3, d4, d5, .false., o1, o2, o3)
  end subroutine primal

  subroutine otis(z, incr, der)
    real(kind=wp), intent(in) :: z(ns, nv)
    real(kind=wp), intent(out) :: incr(ns), der(ns, nv)
    call unpack(z)
    call computeFlowRule_otis(tau, temp, cdeg, back, tssd, tcs, adi, assd, tgnd, agnd, dtime, &
                              incr, rtdi, d1, d2, d3, d4, d5, .true., o1, o2, o3)
    der(:, 1) = d1; der(:, 2) = d2; der(:, 3) = d3; der(:, 4) = d4; der(:, 5) = d5
  end subroutine otis

end program claim2_flowrule
