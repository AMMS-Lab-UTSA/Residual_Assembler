!=======================================================================
!  umat_driver.f90
!
!  Standalone material-point driver for the Oxford / Grilli crystal
!  plasticity UMAT (sources/permissive/ngrilli_Oxford_Crystal_Plasticity).
!
!  Purpose
!  -------
!  Call the Abaqus UMAT OUTSIDE Abaqus for a SINGLE integration point,
!  marching the deformation-gradient history increment-by-increment.
!  Because the UMAT state (STATEV + the /UMPS/ common block) is
!  history-dependent, the whole history is replayed inside ONE process:
!  STATEV and the common block persist from one increment to the next.
!  NEVER jump straight to the final increment.
!
!  It is deliberately agnostic to WHICH `subroutine UMAT(...)` it links
!  against:
!     * aba_stubs/umat_mock.f90  -> a trivial isotropic-elastic UMAT used
!       to exercise the plumbing end-to-end without ifort/Abaqus.
!     * the real umat.for        -> requires ifort + Abaqus (MKL) or a
!       patched source tree; see README.md and build.sh for status.
!
!  I/O format is plain, list-directed (whitespace/newline separated) text
!  so it is trivial to emit/parse from Python (umat_replay.py).  See
!  README.md section "Driver I/O contract" for the byte-for-byte layout.
!
!  Build / run:  see build.sh / build.bat and README.md.
!
!  This driver is original work (Residual_Assembler, MIT).  It compiles
!  against the UMAT sources by path; it does not copy or modify them.
!=======================================================================
program umat_driver
   implicit none

   ! ---- command line ----
   character(len=1024) :: infile, outfile
   integer :: nargs

   ! ---- header ----
   integer :: nstatv, nprops, ntens, ninc
   integer :: ndi, nshr

   ! ---- UMAT argument list (Abaqus standard, real*8 = real(8)) ----
   real(8), allocatable :: stress(:)          ! (ntens)
   real(8), allocatable :: statev(:)          ! (nstatv)
   real(8), allocatable :: ddsdde(:,:)        ! (ntens,ntens)
   real(8), allocatable :: ddsddt(:), drplde(:)
   real(8), allocatable :: stran(:), dstran(:)
   real(8), allocatable :: props(:)           ! (nprops)
   real(8) :: sse, spd, scd, rpl, drpldt, drpldt_dummy
   real(8) :: dtime, temp, dtemp, pnewdt, celent
   real(8) :: time(2), predef(1), dpred(1)
   real(8) :: coords(3), drot(3,3), dfgrd0(3,3), dfgrd1(3,3)
   character(len=80) :: cmname
   integer :: noel, npt, layer, kspt, kstep, kinc

   ! ---- scratch ----
   real(8) :: ftmp(9)
   integer :: inc, i, j, k, uin, uout, ios

   external UMAT

   !--------------------------------------------------------------------
   ! parse args:  umat_driver <input.txt> <output.txt>
   !--------------------------------------------------------------------
   nargs = command_argument_count()
   if (nargs < 2) then
      write(*,'(A)') 'usage: umat_driver <input.txt> <output.txt>'
      write(*,'(A)') '  see README.md for the input/output file format.'
      stop 2
   end if
   call get_command_argument(1, infile)
   call get_command_argument(2, outfile)

   uin = 21
   uout = 22
   open(unit=uin, file=trim(infile), status='old', action='read', iostat=ios)
   if (ios /= 0) then
      write(*,'(A)') 'ERROR: cannot open input file: '//trim(infile)
      stop 3
   end if

   !--------------------------------------------------------------------
   ! header block
   !   line: NSTATV NPROPS NTENS NINC
   !   then: PROPS(1..NPROPS)
   !   then: STATEV0(1..NSTATV)   (initial state; usually all 0 so the
   !                               UMAT kinc<=1 init block fills it)
   !--------------------------------------------------------------------
   read(uin,*,iostat=ios) nstatv, nprops, ntens, ninc
   if (ios /= 0) then
      write(*,'(A)') 'ERROR: failed to read header (NSTATV NPROPS NTENS NINC).'
      stop 4
   end if

   ndi = 3
   nshr = ntens - 3

   allocate(stress(ntens), statev(nstatv), ddsdde(ntens,ntens))
   allocate(ddsddt(ntens), drplde(ntens), stran(ntens), dstran(ntens))
   allocate(props(nprops))

   read(uin,*,iostat=ios) (props(i), i=1,nprops)
   if (ios /= 0) then; write(*,'(A)') 'ERROR: failed to read PROPS.'; stop 5; end if
   read(uin,*,iostat=ios) (statev(i), i=1,nstatv)
   if (ios /= 0) then; write(*,'(A)') 'ERROR: failed to read STATEV0.'; stop 6; end if

   ! constants that Abaqus would supply; harmless for this UMAT
   cmname = 'CPMATERIAL'
   coords = 0.0d0
   drot = 0.0d0
   do i = 1, 3
      drot(i,i) = 1.0d0
   end do
   stran = 0.0d0
   dstran = 0.0d0
   predef = 0.0d0
   dpred = 0.0d0
   celent = 1.0d0
   layer = 1
   kspt = 1
   sse = 0.0d0; spd = 0.0d0; scd = 0.0d0
   rpl = 0.0d0; drpldt = 0.0d0

   open(unit=uout, file=trim(outfile), status='replace', action='write')
   write(uout,'(A)') '# umat_driver output: one record per increment'
   write(uout,'(A)') '# record = KINC / STRESS(ntens) / STATEV(nstatv) / '// &
                     'DDSDDE(ntens*ntens,row-major) / PNEWDT'
   write(uout,'(3(I8,1X))') nstatv, ntens, ninc

   !--------------------------------------------------------------------
   ! increment loop.  Per increment the input provides:
   !   NOEL NPT KSTEP KINC
   !   DTIME TEMP DTEMP TIME1 TIME2
   !   DFGRD0 : 9 values, row-major  (F11 F12 F13 F21 F22 F23 F31 F32 F33)
   !   DFGRD1 : 9 values, row-major
   !--------------------------------------------------------------------
   do inc = 1, ninc

      read(uin,*,iostat=ios) noel, npt, kstep, kinc
      if (ios /= 0) then
         write(*,'(A,I0)') 'ERROR: failed reading inc header at inc=', inc
         stop 7
      end if
      read(uin,*,iostat=ios) dtime, temp, dtemp, time(1), time(2)
      if (ios /= 0) then; write(*,'(A,I0)') 'ERROR: reading inc scalars at inc=', inc; stop 8; end if

      read(uin,*,iostat=ios) (ftmp(k), k=1,9)
      if (ios /= 0) then; write(*,'(A,I0)') 'ERROR: reading DFGRD0 at inc=', inc; stop 9; end if
      call unpack33(ftmp, dfgrd0)

      read(uin,*,iostat=ios) (ftmp(k), k=1,9)
      if (ios /= 0) then; write(*,'(A,I0)') 'ERROR: reading DFGRD1 at inc=', inc; stop 10; end if
      call unpack33(ftmp, dfgrd1)

      pnewdt = 1.0d0
      ddsdde = 0.0d0

      !----------------------------------------------------------------
      ! THE call.  Argument list matches umat.for exactly:
      !   UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
      !        RPL,DDSDDT,DRPLDE,DRPLDT,
      !        STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
      !        NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
      !        CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
      !----------------------------------------------------------------
      call UMAT(stress, statev, ddsdde, sse, spd, scd, &
                rpl, ddsddt, drplde, drpldt, &
                stran, dstran, time, dtime, temp, dtemp, predef, dpred, cmname, &
                ndi, nshr, ntens, nstatv, props, nprops, coords, drot, pnewdt, &
                celent, dfgrd0, dfgrd1, noel, npt, layer, kspt, kstep, kinc)

      ! ---- write this increment's result ----
      write(uout,'(I10)') kinc
      write(uout,'(6(1PE24.16,1X))') (stress(i), i=1,ntens)
      write(uout,'(4(1PE24.16,1X))') (statev(i), i=1,nstatv)
      write(uout,'(6(1PE24.16,1X))') ((ddsdde(i,j), j=1,ntens), i=1,ntens)
      write(uout,'(1PE24.16)') pnewdt

      if (pnewdt < 1.0d0) then
         write(*,'(A,I0,A,1PE12.4)') 'WARNING inc=', inc, &
            ': UMAT requested cutback pnewdt=', pnewdt
      end if
   end do

   close(uin)
   close(uout)
   write(*,'(A,I0,A)') 'umat_driver: wrote ', ninc, ' increment record(s).'

contains

   ! row-major 9-vector -> 3x3 matrix  A(i,j)
   subroutine unpack33(v, A)
      real(8), intent(in)  :: v(9)
      real(8), intent(out) :: A(3,3)
      A(1,1)=v(1); A(1,2)=v(2); A(1,3)=v(3)
      A(2,1)=v(4); A(2,2)=v(5); A(2,3)=v(6)
      A(3,1)=v(7); A(3,2)=v(8); A(3,3)=v(9)
   end subroutine unpack33

end program umat_driver
