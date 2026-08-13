! IIc piece (a): 2-dim double-exponential quadrature of Ga(t,u)
! (g2_iic_a_ga.inc, generated from g2_iic_a_Ga.pkl), double precision.
!
! A* = int_0^1 int_0^1 Ga dt du ~ 0.7697 (Fortran a-column, lam->0).
! Prints A* at two DE levels to show converged digits.
module dli2_mod
   implicit none
   integer, parameter :: dp = selected_real_kind(30)
   real(dp), parameter :: pi = 4*atan(1.0_dp)
   interface cdli2
      module procedure cdli2_c, cdli2_r
   end interface
contains

   complex(dp) function cdli2_r(n, x)
      ! real argument promoted to complex (result may still be complex)
      real(8), intent(in) :: n
      real(dp), intent(in) :: x
      cdli2_r = cdli2_c(n, cmplx(x, 0.0_dp, dp))
   end function

   complex(dp) function cdli2_c(n, z)
      ! complex dilogarithm Li_2(z); n must be 2
      real(8), intent(in) :: n   ! fcode emits a double literal
      complex(dp), intent(in) :: z
      complex(dp) :: w, u, extra, su, up
      real(dp) :: sgn
      integer :: k
      real(dp), save :: b(0:60)
      logical, save :: binit = .false.
      if (.not. binit) then
         call init_b(b)
         binit = .true.
      end if
      if (nint(n) /= 2) error stop "cdli2: n /= 2"
      if (z == (0.0_dp, 0.0_dp)) then
         cdli2_c = 0
         return
      end if
      w = z
      sgn = 1
      extra = 0
      if (abs(w) > 1.0_dp) then
         ! Li2(z) = -Li2(1/z) - pi^2/6 - log(-z)^2/2
         extra = extra - pi**2/6 - 0.5_dp*log(-w)**2
         w = 1/w
         sgn = -sgn
      end if
      if (real(w) > 0.5_dp) then
         ! Li2(w) = pi^2/6 - log(w) log(1-w) - Li2(1-w)
         extra = extra + sgn*(pi**2/6 - log(w)*log(1 - w))
         w = 1 - w
         sgn = -sgn
      end if
      u = -log(1 - w)
      su = 0
      up = u
      do k = 0, 60
         if (b(k) /= 0) su = su + b(k)*up
         up = up*u
      end do
      cdli2_c = sgn*su + extra
   end function

   subroutine init_b(b)
      ! b(k) = B_k/(k+1)! via the exact Bernoulli recursion
      real(dp), intent(out) :: b(0:60)
      real(dp) :: bern(0:60), c, fact
      integer :: m, j
      bern(0) = 1
      do m = 1, 60
         c = 0
         do j = 0, m - 1
            c = c + binom(m + 1, j)*bern(j)
         end do
         bern(m) = -c/(m + 1)
      end do
      fact = 1
      do m = 0, 60
         fact = fact*(m + 1)      ! (m+1)!
         b(m) = bern(m)/fact
      end do
   end subroutine

   real(dp) function binom(n, k)
      integer, intent(in) :: n, k
      integer :: i
      binom = 1
      do i = 1, k
         binom = binom*(n - k + i)/real(i, dp)
      end do
   end function

end module

module gc_mod
   use dli2_mod
   implicit none
   interface logc
      module procedure logc_r, logc_c
   end interface
contains

   complex(dp) function logc_r(x)
      real(dp), intent(in) :: x
      logc_r = log(cmplx(x, 0.0_dp, dp))
   end function

   complex(dp) function logc_c(z)
      complex(dp), intent(in) :: z
      logc_c = log(z)
   end function

   real(dp) function gcfun(s, t, u)
      real(dp), intent(in) :: s, t, u
      real(dp) :: R
      complex(dp) :: gcval
      complex(dp), parameter :: CI = (0.0_dp, 1.0_dp)
      R = sqrt(t*u*(1 - u))
      include "g2_iic_c_gc.inc"
      gcfun = real(gcval, dp)
   end function

end module

program g2_iic_c_quad5
   use gc_mod
   implicit none
   ! C* = int_0^1 dt int_0^{1-t} ds int_0^1 du Gc(s,t,u)
   ! u -> 1 sliver cut at delta, extrapolated (as for piece (a))
   integer, parameter :: ND = 2
   real(dp), parameter :: d0 = 1.6e-5_dp
   real(dp) :: h, res(ND), deltas(ND), A0
   integer :: lev, k

   do lev = 5, 5
      h = 1.0_dp/2**lev
      do k = 1, ND
         deltas(k) = d0*2**(k-1)
         res(k) = de3d(h, deltas(k))
      end do
      ! pure delta^2 Richardson from the first pair; report ladder
      A0 = (4*res(1) - res(2))/3
      print "(a, i2, a, f38.32)", "level ", lev, ": C* (Richardson) = ", A0
      do k = 1, ND
         print "(a, es9.2, a, f38.32)", "   delta=", deltas(k), "  C=", res(k)
      end do
   end do

contains

   real(dp) function de3d(h, delta)
      real(dp), intent(in) :: h, delta
      real(dp), allocatable :: xs(:), ws(:)
      real(dp) :: acc, ti, tj, sc, tv, sv
      integer :: n, i, j, k2
      integer, save :: nbad = 0
      call de_nodes(h, xs, ws, n)
      sc = 1 - delta
      acc = 0
      !$omp parallel do private(i,j,k2,ti,tj,tv,sv) reduction(+:acc) schedule(dynamic)
      do i = 1, n
         tv = xs(i)
         ! Gc evaluation breaks down (cancellation) for t < ~1e-17;
         ! the dropped mass is ~1e-15 (integrand grows only like log^2 t)
         if (tv < 1e-18_dp) cycle
         ti = 0
         do j = 1, n
            sv = (1 - tv)*xs(j)
            tj = 0
            do k2 = 1, n
               block
                  real(dp) :: gv
                  gv = gcfun(sv, tv, sc*xs(k2))
                  if (gv /= gv .or. abs(gv) > 1e8_dp) then
                     !$omp critical
                     nbad = nbad + 1
                     if (nbad <= 10) print "(a,3es12.3,es12.3)", &
                        "bad s,t,u,val: ", sv, tv, sc*xs(k2), gv
                     !$omp end critical
                     gv = 0
                  end if
                  tj = tj + ws(k2)*gv
               end block
            end do
            ti = ti + ws(j)*(1 - tv)*tj
         end do
         acc = acc + ws(i)*ti
      end do
      de3d = acc*sc
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
         if (x < 1e-30_dp .or. 1 - x < 1e-30_dp .or. w < 1e-3000_dp) cycle
         m = m + 1
         xs(m) = x
         ws(m) = w
      end do
      n = m
   end subroutine

end program
