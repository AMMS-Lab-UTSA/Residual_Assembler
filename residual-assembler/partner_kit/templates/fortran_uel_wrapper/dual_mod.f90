! dual_mod.f90 — minimal first-order dual number for Fortran (option 2).
!
! Rewrite your element/global residual against type(dual). Seed a parameter with
! dual(value, 1.0d0); the %du component of the residual is dR/dparam = a column of
! R^(1). This mirrors the C++ Dual1 and the Python Dual1 in the kit.
!
! This is a starter skeleton — extend the overloads your residual actually needs.

module dual_mod
  implicit none
  private
  public :: dual, assignment(=), operator(+), operator(-), operator(*), &
            operator(/), dual_pow, seed, real_part, imag_part

  type :: dual
     real(8) :: re = 0.0d0     ! value
     real(8) :: du = 0.0d0     ! first-order derivative component (eps)
  end type dual

  interface assignment(=)
     module procedure from_real
  end interface

  interface operator(+)
     module procedure add_dd, add_dr, add_rd
  end interface
  interface operator(-)
     module procedure sub_dd, neg_d
  end interface
  interface operator(*)
     module procedure mul_dd, mul_dr, mul_rd
  end interface
  interface operator(/)
     module procedure div_dd
  end interface

contains
  elemental subroutine from_real(a, r)
     type(dual), intent(out) :: a
     real(8), intent(in) :: r
     a%re = r; a%du = 0.0d0
  end subroutine

  elemental type(dual) function seed(v) result(a)
     real(8), intent(in) :: v
     a%re = v; a%du = 1.0d0
  end function

  elemental real(8) function real_part(a) result(r)
     type(dual), intent(in) :: a
     r = a%re
  end function
  elemental real(8) function imag_part(a) result(r)
     type(dual), intent(in) :: a
     r = a%du
  end function

  elemental type(dual) function add_dd(a, b) result(c)
     type(dual), intent(in) :: a, b
     c%re = a%re + b%re; c%du = a%du + b%du
  end function
  elemental type(dual) function add_dr(a, b) result(c)
     type(dual), intent(in) :: a; real(8), intent(in) :: b
     c%re = a%re + b; c%du = a%du
  end function
  elemental type(dual) function add_rd(b, a) result(c)
     real(8), intent(in) :: b; type(dual), intent(in) :: a
     c%re = a%re + b; c%du = a%du
  end function
  elemental type(dual) function sub_dd(a, b) result(c)
     type(dual), intent(in) :: a, b
     c%re = a%re - b%re; c%du = a%du - b%du
  end function
  elemental type(dual) function neg_d(a) result(c)
     type(dual), intent(in) :: a
     c%re = -a%re; c%du = -a%du
  end function
  elemental type(dual) function mul_dd(a, b) result(c)
     type(dual), intent(in) :: a, b
     c%re = a%re * b%re; c%du = a%re * b%du + a%du * b%re
  end function
  elemental type(dual) function mul_dr(a, b) result(c)
     type(dual), intent(in) :: a; real(8), intent(in) :: b
     c%re = a%re * b; c%du = a%du * b
  end function
  elemental type(dual) function mul_rd(b, a) result(c)
     real(8), intent(in) :: b; type(dual), intent(in) :: a
     c%re = a%re * b; c%du = a%du * b
  end function
  elemental type(dual) function div_dd(a, b) result(c)
     type(dual), intent(in) :: a, b
     c%re = a%re / b%re
     c%du = (a%du * b%re - a%re * b%du) / (b%re * b%re)
  end function

  elemental type(dual) function dual_pow(a, n) result(c)
     type(dual), intent(in) :: a; integer, intent(in) :: n
     real(8) :: ap
     ap = a%re ** (n - 1)
     c%re = ap * a%re; c%du = real(n, 8) * ap * a%du
  end function
end module dual_mod
