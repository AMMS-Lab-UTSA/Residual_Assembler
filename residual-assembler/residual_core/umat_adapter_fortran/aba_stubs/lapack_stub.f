c=======================================================================
c  lapack_stub.f  --  portable reference DGETRF / DGETRI
c
c  utils.f::lapinverse() (called by umat.for to invert the 3x3
c  deformation gradient, and by kmat.f) uses LAPACK's DGETRF + DGETRI.
c  Inside Abaqus these come from Intel MKL (the env file compiles with
c  /Qmkl).  Outside Abaqus, if you do NOT link a real LAPACK/MKL, link
c  this file: it supplies mathematically correct (unblocked) DGETRF and
c  DGETRI so lapinverse returns the true inverse.
c
c  These are genuine implementations, NOT fakes:
c    DGETRF -> right-looking LU with partial pivoting (LAPACK IPIV
c              convention: row k is interchanged with row IPIV(k)).
c    DGETRI -> forms inv(A) by solving A*X=I column-by-column with the
c              stored LU factors (equivalent to DGETRS, TRANS='N').
c  They are O(n^3) and unoptimised; for the 3x3 inverses in this UMAT
c  that is irrelevant.  For production accuracy/speed, link real MKL.
c
c  If you DO link real LAPACK/MKL, do NOT also link this file (duplicate
c  symbol).  Residual_Assembler, MIT.
c=======================================================================

      subroutine DGETRF(M, N, A, LDA, IPIV, INFO)
      implicit none
      integer M, N, LDA, INFO
      integer IPIV(*)
      double precision A(LDA,*)
      integer i, j, k, p
      double precision amax, t
      INFO = 0
      do k = 1, min(M,N)
         ! partial pivot: largest |A(i,k)| for i>=k
         p = k
         amax = abs(A(k,k))
         do i = k+1, M
            if (abs(A(i,k)) .gt. amax) then
               amax = abs(A(i,k))
               p = i
            end if
         end do
         IPIV(k) = p
         if (A(p,k) .eq. 0.0d0) then
            if (INFO .eq. 0) INFO = k
         else
            if (p .ne. k) then
               do j = 1, N
                  t = A(k,j)
                  A(k,j) = A(p,j)
                  A(p,j) = t
               end do
            end if
            do i = k+1, M
               A(i,k) = A(i,k) / A(k,k)
            end do
            do j = k+1, N
               do i = k+1, M
                  A(i,j) = A(i,j) - A(i,k)*A(k,j)
               end do
            end do
         end if
      end do
      return
      end

      subroutine DGETRI(N, A, LDA, IPIV, WORK, LWORK, INFO)
      implicit none
      integer N, LDA, LWORK, INFO
      integer IPIV(*)
      double precision A(LDA,*), WORK(*)
      double precision B(N,N), x(N)
      double precision s
      integer i, j, k
      INFO = 0
      if (LWORK .gt. 0) WORK(1) = dble(N)
      ! Solve A * X = I, column by column, using stored LU (A = P L U).
      do k = 1, N
         ! right-hand side = e_k, then apply row pivots (forward)
         do i = 1, N
            x(i) = 0.0d0
         end do
         x(k) = 1.0d0
         do i = 1, N
            call dswap1(x(i), x(IPIV(i)))
         end do
         ! forward solve  L y = Pb   (L unit-lower, stored below diag)
         do i = 1, N
            s = x(i)
            do j = 1, i-1
               s = s - A(i,j)*x(j)
            end do
            x(i) = s
         end do
         ! back solve  U xcol = y     (U on/above diag)
         do i = N, 1, -1
            s = x(i)
            do j = i+1, N
               s = s - A(i,j)*x(j)
            end do
            if (A(i,i) .eq. 0.0d0) then
               INFO = i
               x(i) = 0.0d0
            else
               x(i) = s / A(i,i)
            end if
         end do
         do i = 1, N
            B(i,k) = x(i)
         end do
      end do
      ! copy inverse back into A
      do j = 1, N
         do i = 1, N
            A(i,j) = B(i,j)
         end do
      end do
      return
      end

      subroutine dswap1(a, b)
      implicit none
      double precision a, b, t
      t = a
      a = b
      b = t
      return
      end
