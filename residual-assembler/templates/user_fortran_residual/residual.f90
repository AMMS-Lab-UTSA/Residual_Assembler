! residual.f90 — residual template (Path B: local compiled code).
!
! One binary, two jobs:
!
!   ./residual
!       Self-test. Prints R, the tangent and dR/dk at a known point and
!       cross-checks the analytic parameter derivatives against finite
!       differences of eval_residual. Run it after every edit.
!
!   ./residual --request <req.json> --response <resp.npz>
!       Black-box responder. This is what `resasm run` invokes (see
!       resasm.yml -> residual.command). It reads the framework's request,
!       computes the tangent dR/du and the order-1 right-hand side dR/da_i,
!       and writes the response file.
!
! RESPONSE FILE NOTE: the framework hands us a --response path ending in `.npz`.
! A compiled program has no business writing numpy archives, so we strip the
! trailing `.npz` and write `<that>.json` instead. The framework's reader looks
! for exactly that fallback (see resasm_user/providers.py::_read_response).
!
! For arbitrary-order sensitivities, rewrite the model against OTILib's Fortran
! OTI type (F95+, static-dense) — see ../../partner_kit/templates/.

! ===========================================================================
! 1. YOUR MODEL — this is the only module you need to edit.
!
!    R = k u^3 - f        (1-DOF cubic spring)
!    params(1) = k, params(2) = f   (same order as `parameters:` in resasm.yml)
! ===========================================================================
module residual_mod
  implicit none
  integer, parameter :: MODEL_NDOF   = 1   ! unknowns      (problem.unknowns)
  integer, parameter :: MODEL_NPARAM = 2   ! parameters, in resasm.yml key order
contains

  ! The residual itself.
  subroutine eval_residual(u, params, ndof, R)
    integer, intent(in) :: ndof
    double precision, intent(in)  :: u(ndof), params(MODEL_NPARAM)
    double precision, intent(out) :: R(ndof)
    R(1) = params(1) * u(1)**3 - params(2)
  end subroutine eval_residual

  ! The tangent dR/du. Used for the sensitivity solve.
  subroutine eval_tangent(u, params, ndof, T)
    integer, intent(in) :: ndof
    double precision, intent(in)  :: u(ndof), params(MODEL_NPARAM)
    double precision, intent(out) :: T(ndof, ndof)
    T(1,1) = 3.0d0 * params(1) * u(1)**2
  end subroutine eval_tangent

  ! dR/d(parameter ip), analytic. `ip` is 1-based, in resasm.yml parameter order.
  !   dR/dk = u^3        dR/df = -1
  ! EDIT THIS WHENEVER YOU EDIT eval_residual. The self-test below finite-
  ! differences eval_residual and will fail loudly if the two ever drift apart.
  subroutine eval_dR_dparam(u, params, ndof, ip, dR)
    integer, intent(in) :: ndof, ip
    double precision, intent(in)  :: u(ndof), params(MODEL_NPARAM)
    double precision, intent(out) :: dR(ndof)
    if (ip == 1) then
      dR(1) = u(1)**3                 ! dR/dk
    else
      dR(1) = -1.0d0                  ! dR/df
    end if
  end subroutine eval_dR_dparam

end module residual_mod


! ===========================================================================
! 2. Minimal JSON reading.
!
! The request has a fixed, simple shape, so a tolerant scanner is plenty: we
! only need "u":[...], "parameters":{...}, "order":N and
! "direction_map":{"1":[[...],[...]]}. No escapes, no nesting surprises.
! ===========================================================================
module json_mod
  implicit none
  integer, parameter :: MAXN = 64          ! cap on parsed vector / column counts
contains

  ! Read a whole file into a string.
  subroutine slurp(path, buf)
    character(len=*), intent(in) :: path
    character(len=:), allocatable, intent(out) :: buf
    integer :: iu, sz, ios
    logical :: ex
    inquire(file=path, exist=ex, size=sz)
    if (.not. ex .or. sz <= 0) then
      write(0,'(a,a)') 'residual: cannot open request file: ', trim(path)
      stop 2
    end if
    allocate(character(len=sz) :: buf)
    open(newunit=iu, file=path, access='stream', form='unformatted', &
         status='old', action='read', iostat=ios)
    if (ios /= 0) then
      write(0,'(a,a)') 'residual: cannot open request file: ', trim(path)
      stop 2
    end if
    read(iu) buf
    close(iu)
  end subroutine slurp

  ! Index of the first character of the value of "key", searching from `from`.
  ! Returns 0 if the key is not present.
  integer function value_pos(s, key, from)
    character(len=*), intent(in) :: s, key
    integer, intent(in) :: from
    integer :: p, c
    value_pos = 0
    p = index(s(from:), '"'//key//'"')
    if (p == 0) return
    p = from + p - 1 + len(key) + 2        ! first char after the closing quote
    c = index(s(p:), ':')
    if (c == 0) return
    value_pos = p + c                      ! first char after the ':'
  end function value_pos

  ! Skip whitespace and the separators we never care about.
  subroutine skip_sep(s, p)
    character(len=*), intent(in) :: s
    integer, intent(inout) :: p
    do while (p <= len(s))
      if (index(' '//achar(9)//achar(10)//achar(13)//',', s(p:p)) == 0) exit
      p = p + 1
    end do
  end subroutine skip_sep

  ! One number at/after p; p advances past it.
  double precision function read_number(s, p)
    character(len=*), intent(in) :: s
    integer, intent(inout) :: p
    integer :: i, j, n
    n = len(s)
    read_number = 0.0d0
    call skip_sep(s, p)
    i = p
    j = i
    do while (j <= n)
      if (index('+-0123456789.eE', s(j:j)) == 0) exit
      j = j + 1
    end do
    if (j > i) read(s(i:j-1), *) read_number
    p = j
  end function read_number

  ! "[a, b, c]" at/after p -> vals(1..n); p ends just past the closing ']'.
  subroutine read_number_array(s, p, vals, n)
    character(len=*), intent(in) :: s
    integer, intent(inout) :: p
    double precision, intent(out) :: vals(:)
    integer, intent(out) :: n
    n = 0
    do while (p <= len(s))
      if (s(p:p) == '[') exit
      p = p + 1
    end do
    p = p + 1                              ! past '['
    do
      call skip_sep(s, p)
      if (p > len(s)) exit
      if (s(p:p) == ']') then
        p = p + 1
        exit
      end if
      if (n >= size(vals)) then
        write(0,'(a)') 'residual: JSON array longer than MAXN'
        stop 2
      end if
      n = n + 1
      vals(n) = read_number(s, p)
    end do
  end subroutine read_number_array

  ! {"k": 2.0, "f": 16.0} -> names + values IN FILE ORDER, which is the
  ! resasm.yml parameter order (the request's "seed_directions" says the same).
  subroutine read_number_object(s, p, names, vals, n)
    character(len=*), intent(in) :: s
    integer, intent(inout) :: p
    character(len=*), intent(out) :: names(:)
    double precision, intent(out) :: vals(:)
    integer, intent(out) :: n
    integer :: q, c
    n = 0
    do while (p <= len(s))
      if (s(p:p) == '{') exit
      p = p + 1
    end do
    p = p + 1                              ! past '{'
    do
      call skip_sep(s, p)
      if (p > len(s)) exit
      if (s(p:p) == '}') then
        p = p + 1
        exit
      end if
      if (s(p:p) /= '"') then
        p = p + 1
        cycle
      end if
      q = p + index(s(p+1:), '"')          ! closing quote of the name
      if (n >= size(vals)) then
        write(0,'(a)') 'residual: JSON object larger than MAXN'
        stop 2
      end if
      n = n + 1
      names(n) = s(p+1:q-1)
      c = q + index(s(q:), ':')             ! first char after the ':'
      p = c
      vals(n) = read_number(s, p)
    end do
  end subroutine read_number_object

  ! "direction_map": {"1": [[1,0],[0,1]]} -> the exponent vectors for `order`,
  ! one per column of the order-`order` right-hand side.
  subroutine read_direction_map(s, order, exps, m, ncols)
    character(len=*), intent(in) :: s
    integer, intent(in) :: order
    integer, intent(out) :: exps(:,:)      ! (m, ncols)
    integer, intent(out) :: m, ncols
    integer :: dm, p, n, i
    character(len=16) :: okey
    double precision :: tmp(size(exps,1))
    ncols = 0
    m = 0
    dm = value_pos(s, 'direction_map', 1)
    if (dm == 0) return
    write(okey, '(i0)') order
    p = value_pos(s, trim(okey), dm)       ! the "<order>" key inside the map
    if (p == 0) return
    do while (p <= len(s))
      if (s(p:p) == '[') exit
      p = p + 1
    end do
    p = p + 1                              ! past the OUTER '['
    do
      call skip_sep(s, p)
      if (p > len(s)) exit
      if (s(p:p) == ']') then
        p = p + 1
        exit
      end if
      call read_number_array(s, p, tmp, n) ! one exponent vector
      if (ncols >= size(exps, 2)) then
        write(0,'(a)') 'residual: more directions than MAXN'
        stop 2
      end if
      ncols = ncols + 1
      m = n
      do i = 1, n
        exps(i, ncols) = nint(tmp(i))
      end do
    end do
  end subroutine read_direction_map

end module json_mod


! ===========================================================================
! 3. Driver: self-test (no args) or black-box responder (--request/--response).
! ===========================================================================
program main
  use residual_mod
  use json_mod
  implicit none

  character(len=4096) :: arg, reqf, respf
  integer :: i, nargs, st
  logical :: hasreq, hasresp

  hasreq  = .false.
  hasresp = .false.
  reqf    = ''
  respf   = ''
  nargs   = command_argument_count()

  i = 1
  do while (i < nargs)
    call get_command_argument(i, arg)
    if (trim(arg) == '--request') then
      call get_command_argument(i + 1, reqf)
      hasreq = .true.
      i = i + 2
    else if (trim(arg) == '--response') then
      call get_command_argument(i + 1, respf)
      hasresp = .true.
      i = i + 2
    else
      i = i + 1
    end if
  end do

  if (hasreq .and. hasresp) then
    call respond(trim(reqf), trim(respf), st)      ! black-box responder
    if (st /= 0) stop 2
  else if (hasreq .or. hasresp) then
    write(0,'(a)') 'usage: residual [--request <in.json> --response <out.npz>]'
    stop 2
  else
    call self_test(st)                             ! no args -> self-test
    if (st /= 0) stop 1
  end if

contains

  ! ---- black-box responder ------------------------------------------------
  subroutine respond(reqf, respf, st)
    character(len=*), intent(in) :: reqf, respf
    integer, intent(out) :: st

    character(len=:), allocatable :: s, outp
    double precision :: u(MAXN), pv(MAXN)
    character(len=64) :: pnames(MAXN)
    integer :: exps(MAXN, MAXN)
    double precision :: T(MODEL_NDOF, MODEL_NDOF), col(MODEL_NDOF)
    double precision, allocatable :: R1(:,:)
    integer :: p, nu, npar, order, m, ncols, c, ip, e, i

    st = 0
    call slurp(reqf, s)

    ! --- read the request --------------------------------------------------
    p = value_pos(s, 'u', 1)
    call read_number_array(s, p, u, nu)

    p = value_pos(s, 'parameters', 1)
    call read_number_object(s, p, pnames, pv, npar)

    p = value_pos(s, 'order', 1)
    order = nint(read_number(s, p))

    call read_direction_map(s, order, exps, m, ncols)

    ! --- sanity: the model routines above are written for one specific shape.
    !     Say so rather than silently returning nonsense. ---------------------
    if (nu /= MODEL_NDOF .or. npar /= MODEL_NPARAM) then
      write(0,'(a,i0,a,i0,a,i0,a,i0,a)') &
        'residual: request has ndof=', nu, ', nparam=', npar, &
        ' but this model is written for ndof=', MODEL_NDOF, &
        ', nparam=', MODEL_NPARAM, '. Update MODEL_NDOF / MODEL_NPARAM.'
      st = 2
      return
    end if
    if (order /= 1) then
      write(0,'(a,i0,a)') &
        'residual: this template implements order 1 only (got order=', order, &
        '). For order >= 2 you must also consume "u_star_coefficients".'
      st = 2
      return
    end if
    if (ncols == 0) then
      write(0,'(a)') 'residual: request has no direction_map entry for this order.'
      st = 2
      return
    end if

    ! --- compute -----------------------------------------------------------
    call eval_tangent(u(1:MODEL_NDOF), pv(1:MODEL_NPARAM), MODEL_NDOF, T)

    ! order-1 RHS: one column per direction; R1(:,c) = dR/da_ip, where ip is the
    ! position of the 1 in that direction's exponent vector.
    allocate(R1(MODEL_NDOF, ncols))
    do c = 1, ncols
      ip = 1
      do e = 1, m
        if (exps(e, c) /= 0) then
          ip = e
          exit
        end if
      end do
      call eval_dR_dparam(u(1:MODEL_NDOF), pv(1:MODEL_NPARAM), MODEL_NDOF, ip, col)
      do i = 1, MODEL_NDOF
        R1(i, c) = col(i)
      end do
    end do

    ! --- write the response ------------------------------------------------
    outp = response_json_path(respf)
    call write_response(outp, order, R1, MODEL_NDOF, ncols, T)
  end subroutine respond

  ! The framework gives us "<dir>/response.npz"; we write "<dir>/response.json".
  function response_json_path(resp) result(outp)
    character(len=*), intent(in) :: resp
    character(len=:), allocatable :: outp
    integer :: n
    outp = resp
    n = len(outp)
    if (n > 4) then
      if (outp(n-3:n) == '.npz') outp = outp(1:n-4)
    end if
    n = len(outp)
    if (n < 5) then
      outp = outp//'.json'
    else if (outp(n-4:n) /= '.json') then
      outp = outp//'.json'
    end if
  end function response_json_path

  subroutine write_response(path, order, R1, ndof, ncols, T)
    character(len=*), intent(in) :: path
    integer, intent(in) :: order, ndof, ncols
    double precision, intent(in) :: R1(:,:), T(:,:)
    integer :: iu
    character(len=16) :: okey
    open(newunit=iu, file=path, status='replace', action='write')
    write(okey,'(i0)') order
    write(iu,'(a)', advance='no') '{"residual_coefficients_by_order": {"'// &
                                  trim(okey)//'": '
    call write_matrix(iu, R1, ndof, ncols)
    write(iu,'(a)', advance='no') '}, "tangent": '
    call write_matrix(iu, T, ndof, ndof)
    write(iu,'(a)') ', "diagnostics": {"solver": "template_fortran", '// &
                    '"method": "analytic"}}'
    close(iu)
  end subroutine write_response

  subroutine write_matrix(iu, a, nr, nc)
    integer, intent(in) :: iu, nr, nc
    double precision, intent(in) :: a(:,:)
    integer :: i, j
    character(len=32) :: num
    write(iu,'(a)', advance='no') '['
    do i = 1, nr
      if (i > 1) write(iu,'(a)', advance='no') ','
      write(iu,'(a)', advance='no') '['
      do j = 1, nc
        if (j > 1) write(iu,'(a)', advance='no') ','
        write(num,'(es24.16e3)') a(i,j)
        write(iu,'(a)', advance='no') trim(adjustl(num))
      end do
      write(iu,'(a)', advance='no') ']'
    end do
    write(iu,'(a)', advance='no') ']'
  end subroutine write_matrix

  ! ---- self-test ----------------------------------------------------------
  subroutine self_test(st)
    integer, intent(out) :: st
    double precision :: u(MODEL_NDOF), p(MODEL_NPARAM), R(MODEL_NDOF)
    double precision :: T(MODEL_NDOF, MODEL_NDOF), dk(MODEL_NDOF)
    double precision :: pp(MODEL_NPARAM), pm(MODEL_NPARAM)
    double precision :: Rp(MODEL_NDOF), Rm(MODEL_NDOF), an(MODEL_NDOF), h, fd
    integer :: ip, i

    st = 0
    u = 2.0d0
    p = (/ 2.0d0, 16.0d0 /)                        ! k, f
    call eval_residual(u, p, MODEL_NDOF, R)
    call eval_tangent(u, p, MODEL_NDOF, T)
    print *, "R =", R(1), " (expect 0)"
    print *, "T = dR/du =", T(1,1), " (expect 24)"

    call eval_dR_dparam(u, p, MODEL_NDOF, 1, dk)
    print *, "dR/dk =", dk(1), " du/dk =", -dk(1) / T(1,1)

    ! Guard against the analytic dR/da drifting away from eval_residual:
    ! central finite differences must reproduce eval_dR_dparam.
    do ip = 1, MODEL_NPARAM
      h = 1.0d-6 * max(1.0d0, abs(p(ip)))
      pp = p
      pm = p
      pp(ip) = pp(ip) + h
      pm(ip) = pm(ip) - h
      call eval_residual(u, pp, MODEL_NDOF, Rp)
      call eval_residual(u, pm, MODEL_NDOF, Rm)
      call eval_dR_dparam(u, p, MODEL_NDOF, ip, an)
      do i = 1, MODEL_NDOF
        fd = (Rp(i) - Rm(i)) / (2.0d0 * h)
        if (abs(fd - an(i)) > 1.0d-6 * max(1.0d0, abs(fd))) then
          print *, "MISMATCH dR/da: parameter", ip, " analytic", an(i), " FD", fd
          st = 1
        end if
      end do
    end do
    if (st == 0) then
      print *, "analytic dR/da vs finite differences: ok"
    else
      print *, "analytic dR/da vs finite differences: FAIL"
    end if
  end subroutine self_test

end program main
