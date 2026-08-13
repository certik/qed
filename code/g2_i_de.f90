! Diagram I (crossed ladder) at lam=0: 5-dim tanh-sinh product rule.
! mu_I = int fI du dv dr dy dt  (u+v+r<1, y+t<1)
! GL failed here: the integrand has (integrable) boundary singularities
! at v->1 etc.; DE clusters nodes at the boundary.  Double-precision
! evaluation breaks down within ~1e-5 of the simplex faces, so each unit
! coordinate is interval-scaled to (delta, 1-delta) and the missing
! boundary mass is removed by a delta-ladder fit.
! target: 1/6 + 13/36 pi^2 + 5/4 zeta3 - 5/6 pi^2 log2 = -0.4676445...
module i_de_mod
   implicit none
   integer, parameter :: dp = kind(1.d0)
   real(dp), parameter :: pi = 4*atan(1.0_dp)
   integer :: nbad = 0
   interface logc
      module procedure logc_r, logc_c
   end interface
contains
   complex(dp) function logc_r(x)
      real(dp), intent(in) :: x
      logc_r = log(cmplx(x, 0.0_dp, dp))
   end function
   complex(dp) function logc_c(zz)
      complex(dp), intent(in) :: zz
      logc_c = log(zz)
   end function
   real(dp) function ffI(u, v, r, y, t, lam)
      real(dp), intent(in) :: u, v, r, y, t, lam
      complex(dp) :: fIv
      complex(dp), parameter :: CI = (0.0_dp, 1.0_dp)
      ! g2_i_f_split.inc = g2_i_f.inc split into per-term statements by
      ! code/g2_i_split.py (flang chokes on the single 850 KB expression:
      ! >20 min, >30 GB at any -O level; the split compiles in seconds)
      include "g2_i_f_split_decl.inc"
      include "g2_i_f_split.inc"
      ffI = real(fIv, dp)
      if (ffI /= ffI) then
         !$omp atomic
         nbad = nbad + 1
         ffI = 0
      end if
   end function
end module

program g2_i_de
   use i_de_mod
   implicit none
   real(dp), parameter :: deltas(3) = [1e-4_dp, 1e-5_dp, 1e-6_dp]
   real(dp) :: val
   integer :: id, lev
   print "(a)", "  delta      mu_I(delta)        target = -0.467645446094"
   do lev = 2, 3
      print "(a,i2)", "level ", lev
      do id = 1, size(deltas)
         nbad = 0
         val = de5d(1.0_dp/2**lev, deltas(id))
         print "(es9.2, f20.12, a, i9)", deltas(id), val, "   nbad=", nbad
      end do
   end do
contains
   real(dp) function de5d(h, delta)
      real(dp), intent(in) :: h, delta
      real(dp), allocatable :: xs(:), ws(:)
      real(dp) :: acc, a1, a2, a3, a4
      real(dp) :: u, v, r, y, t, jac
      integer :: n, i1, i2, i3, i4, i5
      call de_nodes(h, delta, xs, ws, n)
      acc = 0
      !$omp parallel do private(i1,i2,i3,i4,i5,u,v,r,y,t,jac,a1,a2,a3,a4) &
      !$omp reduction(+:acc) schedule(dynamic)
      do i1 = 1, n
         v = xs(i1)
         a1 = 0
         do i2 = 1, n
            r = (1 - v)*xs(i2)
            a2 = 0
            do i3 = 1, n
               u = (1 - v - r)*xs(i3)
               a3 = 0
               do i4 = 1, n
                  t = xs(i4)
                  a4 = 0
                  do i5 = 1, n
                     y = (1 - t)*xs(i5)
                     a4 = a4 + ws(i5)*ffI(u, v, r, y, t, 0.0_dp)
                  end do
                  a3 = a3 + ws(i4)*(1 - t)*a4
               end do
               a2 = a2 + ws(i3)*a3
            end do
            a1 = a1 + ws(i2)*(1 - v - r)*a2
         end do
         acc = acc + ws(i1)*(1 - v)*a1
      end do
      de5d = acc
   end function
   subroutine de_nodes(h, delta, xs, ws, n)
      ! tanh-sinh nodes on (0,1), affinely squeezed into (delta, 1-delta)
      real(dp), intent(in) :: h, delta
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
         if (x < 1e-17_dp .or. 1 - x < 1e-17_dp .or. w < 1e-300_dp) cycle
         m = m + 1
         xs(m) = delta + (1 - 2*delta)*x
         ws(m) = (1 - 2*delta)*w
      end do
      n = m
   end subroutine
end program
