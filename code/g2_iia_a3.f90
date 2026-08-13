! IIa piece (a): 3-dim quad DE integral of s*chi*g over
! (s<1-t, chi<1, t<1) at a lambda ladder; result - (1/2)log(lam) -> A0
! predicted A0 = 11/12 + pi^2/18 = 1.4649780222837...
module iia_a3_mod
   implicit none
   integer, parameter :: dp = selected_real_kind(30)
   real(dp), parameter :: pi = 4*atan(1.0_dp)
contains
   real(dp) function gfun(s, chi, t, lam)
      real(dp), intent(in) :: s, chi, t, lam
      real(dp) :: gv
      include "g2_iia_a3.inc"
      gfun = gv
   end function
end module

program g2_iia_a3
   use iia_a3_mod
   implicit none
   real(dp), parameter :: lams(6) = [1e-2_dp, 3e-3_dp, 1e-3_dp, 3e-4_dp, 1e-4_dp, 3e-5_dp]
   real(dp) :: lam, val
   integer :: il, lev
   do lev = 6, 7
      print "(a,i2)", "level ", lev
      do il = 1, size(lams)
         lam = lams(il)
         val = de3d(1.0_dp/2**lev, lam)
         print "(es9.2, f36.30)", lam, val - 0.5_dp*log(lam)
      end do
   end do
contains
   real(dp) function de3d(h, lam)
      real(dp), intent(in) :: h, lam
      real(dp), allocatable :: xs(:), ws(:)
      real(dp) :: acc, ti, tj, tv, sv, cv
      integer :: n, i, j, k
      call de_nodes(h, xs, ws, n)
      acc = 0
      !$omp parallel do private(i,j,k,ti,tj,tv,sv,cv) reduction(+:acc) schedule(dynamic)
      do i = 1, n
         tv = xs(i)
         ti = 0
         do j = 1, n
            sv = (1 - tv)*xs(j)
            tj = 0
            do k = 1, n
               cv = xs(k)
               tj = tj + ws(k)*gfun(sv, cv, tv, lam)
            end do
            ti = ti + ws(j)*(1 - tv)*tj
         end do
         acc = acc + ws(i)*ti
      end do
      de3d = acc
   end function
   subroutine de_nodes(h, xs, ws, n)
      real(dp), intent(in) :: h
      real(dp), allocatable, intent(out) :: xs(:), ws(:)
      integer, intent(out) :: n
      real(dp) :: tk, x, w, sh
      integer :: k, kmax, m
      kmax = int(6.0_dp/h)
      allocate(xs(2*kmax + 1), ws(2*kmax + 1))
      m = 0
      do k = -kmax, kmax
         tk = k*h
         sh = 0.5_dp*pi*sinh(tk)
         x = 0.5_dp*(1 + tanh(sh))
         w = 0.25_dp*pi*h*cosh(tk)/cosh(sh)**2
         if (x < 1e-22_dp .or. 1 - x < 1e-22_dp .or. w < 1e-3000_dp) cycle
         m = m + 1
         xs(m) = x
         ws(m) = w
      end do
      n = m
   end subroutine
end program
