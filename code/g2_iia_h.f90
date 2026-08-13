! IIa piece (a): 2-dim quad-precision integral of h(s,t;lam) with a
! lambda ladder; I(lam) - (1/2) log(lam) must converge to A0.
module iia_h_mod
   implicit none
   integer, parameter :: dp = selected_real_kind(30)
   real(dp), parameter :: pi = 4*atan(1.0_dp)
   interface logc
      module procedure logc_r, logc_c
   end interface
   interface atanc
      module procedure atanc_r, atanc_c
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
   complex(dp) function atanc_r(x)
      real(dp), intent(in) :: x
      atanc_r = atan(x)
   end function
   complex(dp) function atanc_c(zz)
      ! complex atan via logs
      complex(dp), intent(in) :: zz
      complex(dp), parameter :: CIu = (0.0_dp, 1.0_dp)
      atanc_c = (CIu/2)*(log(1 - CIu*zz) - log(1 + CIu*zz))
   end function
   real(dp) function hfun(s, t, lam)
      real(dp), intent(in) :: s, t, lam
      complex(dp) :: hv
      complex(dp), parameter :: CI = (0.0_dp, 1.0_dp)
      include "g2_iia_h.inc"
      hfun = real(hv, dp)
   end function
end module

program g2_iia_h
   use iia_h_mod
   implicit none
   real(dp), parameter :: lams(4) = [1e-2_dp, 5e-3_dp, 2e-3_dp, 1e-3_dp]
   real(dp) :: lam, res(3), R1
   integer :: il, lev, k
   do lev = 6, 6
      print "(a,i2)", "level ", lev
      do il = 1, size(lams)
         lam = lams(il)
         do k = 1, 2
            res(k) = de2d(1.0_dp/2**lev, lam, 1.0e-4_dp*2**(k-1))
         end do
         R1 = (4*res(1) - res(2))/3
         print "(es9.2, f36.30)", lam, R1 - 0.5_dp*log(lam)
      end do
   end do
contains
   real(dp) function de2d(h, lam, delta)
      real(dp), intent(in) :: h, lam, delta
      real(dp), allocatable :: xs(:), ws(:)
      real(dp) :: acc, ti, tv, sv
      integer :: n, i, j
      call de_nodes(h, xs, ws, n)
      acc = 0
      !$omp parallel do private(i,j,ti,tv,sv) reduction(+:acc) schedule(dynamic)
      do i = 1, n
         tv = (1 - delta)*xs(i)
         ti = 0
         do j = 1, n
            sv = (1 - tv)*xs(j)
            if (min(sv, tv, 1 - tv, 1 - tv - sv) < 1e-14_dp) cycle
            ti = ti + ws(j)*(1 - tv)*hfun(sv, tv, lam)
         end do
         acc = acc + ws(i)*ti
      end do
      de2d = acc*(1 - delta)
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
         if (x < 1e-25_dp .or. 1 - x < 1e-25_dp .or. w < 1e-3000_dp) cycle
         m = m + 1
         xs(m) = x
         ws(m) = w
      end do
      n = m
   end subroutine
end program
