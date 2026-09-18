! Timing driver for the slide-32 cost comparison (claim 3).
!
! Linked against ONE compiled provider object (umat-oti-provider build), which
! bundles the author's ORIGINAL UMAT and its OTI lift compiled with the same
! flags. Times, for the same strain path at NIP integration points:
!   plain  : the ORIGINAL UMAT, increment by increment        (no sensitivities)
!   eval   : UMAT_OTI_EVAL with the incoming derivatives carried (all parameters)
!   march  : UMAT_OTI_MARCH, the whole path in one call         (all parameters)
! Input on stdin: reps nip npath nprops nparam dtime / props(nprops) / dstran(6).
! Output: one line per method with the wall-clock seconds of ONE analysis
! (reps repetitions after one untimed warm-up repetition).
program cp_timing
  implicit none
  integer :: reps, nip, npath, nprops, nparam, r, ip, n, i
  integer, parameter :: nt = 6, ns = 1
  real(8) :: dtime
  real(8), allocatable :: props(:), path(:,:), dtarr(:), dsig(:,:,:), ddout(:,:,:)
  real(8), allocatable :: dsp(:,:), dstp(:,:), dsp_in(:,:), dstp_in(:,:)
  real(8) :: dstran(nt), stress(nt), statev(ns), ddsdde(nt,nt), stran(nt), strout(nt)
  real(8) :: sse, spd, scd, rpl, ddsddt(nt), drplde(nt), drpldt, time(2), predef(1), dpred(1)
  real(8) :: coords(3), drot(3,3), dfgrd0(3,3), dfgrd1(3,3), pnewdt, celent, checksum
  character(80) :: cmname
  integer(8) :: c0, c1, rate

  read(*,*) reps, nip, npath, nprops, nparam, dtime
  allocate(props(nprops), path(nt,npath), dtarr(npath), dsig(nt,nparam,npath), ddout(nt,nt,npath))
  allocate(dsp(nt,nparam), dstp(ns,nparam), dsp_in(nt,nparam), dstp_in(ns,nparam))
  read(*,*) props
  read(*,*) dstran
  do n = 1, npath
    path(:,n) = dstran
    dtarr(n) = dtime
  end do
  cmname = 'MATERIAL'
  drot = 0.d0; dfgrd0 = 0.d0; dfgrd1 = 0.d0; coords = 0.d0
  do i = 1, 3
    drot(i,i) = 1.d0; dfgrd0(i,i) = 1.d0; dfgrd1(i,i) = 1.d0
  end do
  checksum = 0.d0

  ! ---- plain ORIGINAL UMAT --------------------------------------------------
  do r = 0, reps
    if (r == 1) call system_clock(c0, rate)
    do ip = 1, nip
      stress = 0.d0; statev = 0.d0; stran = 0.d0; time = 0.d0
      do n = 1, npath
        sse = 0.d0; spd = 0.d0; scd = 0.d0; rpl = 0.d0; ddsddt = 0.d0; drplde = 0.d0
        drpldt = 0.d0; pnewdt = 1.d0; celent = 1.d0
        call umat(stress, statev, ddsdde, sse, spd, scd, rpl, ddsddt, drplde, drpldt, &
                  stran, path(:,n), time, dtime, 0.d0, 0.d0, predef, dpred, cmname, &
                  3, 3, nt, ns, props, nprops, coords, drot, pnewdt, celent, dfgrd0, dfgrd1, &
                  1, ip, 1, 1, 1, n)
        stran = stran + path(:,n); time = time + dtime
      end do
      checksum = checksum + stress(4)
    end do
  end do
  call system_clock(c1)
  write(*,'(A,ES24.16)') 'plain ', dble(c1 - c0) / dble(rate) / dble(reps)

  ! ---- OTI EVAL per increment, derivatives carried --------------------------
  do r = 0, reps
    if (r == 1) call system_clock(c0, rate)
    do ip = 1, nip
      stress = 0.d0; statev = 0.d0; stran = 0.d0; time = 0.d0
      dsp_in = 0.d0; dstp_in = 0.d0
      do n = 1, npath
        call umat_oti_eval(stress, statev, ddsdde, stran, path(:,n), time, dtime, 0.d0, 0.d0, &
                           props, nprops, nt, ns, nparam, dsp, dstp, dsp_in, dstp_in)
        dsp_in = dsp; dstp_in = dstp
        stran = stran + path(:,n); time = time + dtime
      end do
      checksum = checksum + dsp(4,1)
    end do
  end do
  call system_clock(c1)
  write(*,'(A,ES24.16)') 'eval  ', dble(c1 - c0) / dble(rate) / dble(reps)

  ! ---- OTI MARCH, whole path in one call ------------------------------------
  do r = 0, reps
    if (r == 1) call system_clock(c0, rate)
    do ip = 1, nip
      call umat_oti_march(props, nprops, path, npath, dtarr, nt, ns, nparam, dsig, strout, ddout)
      checksum = checksum + dsig(4,1,npath)
    end do
  end do
  call system_clock(c1)
  write(*,'(A,ES24.16)') 'march ', dble(c1 - c0) / dble(rate) / dble(reps)
  write(*,'(A,ES24.16)') 'check ', checksum
end program cp_timing
