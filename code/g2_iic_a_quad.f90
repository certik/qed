! IIc piece (a): 2-dim double-exponential quadrature of Ga(t,u)
! (g2_iic_a_ga.inc, generated from g2_iic_a_Ga.pkl), double precision.
!
! A* = int_0^1 int_0^1 Ga dt du ~ 0.7697 (Fortran a-column, lam->0).
! Prints A* at two DE levels to show converged digits.
module dli2_mod
   implicit none
   integer, parameter :: dp = selected_real_kind(30)
   real(dp), parameter :: pi = 4*atan(1.0_dp)
contains

   complex(dp) function cdli2(n, z)
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
         cdli2 = 0
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
      cdli2 = sgn*su + extra
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

module ga_mod
   use dli2_mod
   implicit none
   interface logc
      module procedure logc_r, logc_c
   end interface
contains

   complex(dp) function logc_r(x)
      ! complex-promoting log: log(negative real) = log|x| + i pi
      real(dp), intent(in) :: x
      logc_r = log(cmplx(x, 0.0_dp, dp))
   end function

   complex(dp) function logc_c(z)
      complex(dp), intent(in) :: z
      logc_c = log(z)
   end function

   real(dp) function gafun(t, u)
      real(dp), intent(in) :: t, u
      real(dp) :: R
      complex(dp) :: gaval
      complex(dp), parameter :: CI = (0.0_dp, 1.0_dp)
      R = sqrt(t*u*(1 - u))
      include "g2_iic_a_ga.inc"
      gafun = real(gaval, dp)
   end function

end module

program g2_iic_a_quad
   use ga_mod
   implicit none
   ! Ga has catastrophic cancellation for 1-u < ~1e-4 (double precision);
   ! integrate u over [0, 1-delta] for several delta and extrapolate
   ! delta -> 0 with A(delta) = A0 + a d^2 + b d^2 log d + c d^3 + e d^3 log d
   integer, parameter :: ND = 5
   real(dp), parameter :: d0 = 8.0e-6_dp
   real(dp) :: h, res(ND), deltas(ND), m(5,5), rhs(5), A0
   integer :: lev, k

   do lev = 8, 8
      h = 1.0_dp/2**lev
      do k = 1, ND
         deltas(k) = d0*2**(k-1)
         res(k) = de2d(h, deltas(k))
      end do
      do k = 1, ND
         m(k,:) = [1.0_dp, deltas(k)**2, deltas(k)**2*log(deltas(k)), &
                   deltas(k)**3, deltas(k)**3*log(deltas(k))]
         rhs(k) = res(k)
      end do
      call solve5(m, rhs)
      A0 = rhs(1)
      print "(a, i2, a, f38.32)", "level ", lev, ": A* (extrap) = ", A0
      do k = 1, ND
         print "(a, es9.2, a, f38.32)", "   delta=", deltas(k), "  A=", res(k)
      end do
   end do

contains

   subroutine solve5(a, b)
      real(dp), intent(inout) :: a(5,5), b(5)
      integer :: i, j, p
      real(dp) :: f
      do i = 1, 5
         p = maxloc(abs(a(i:5,i)), 1) + i - 1
         if (p /= i) then
            a([i,p],:) = a([p,i],:)
            b([i,p]) = b([p,i])
         end if
         do j = i+1, 5
            f = a(j,i)/a(i,i)
            a(j,:) = a(j,:) - f*a(i,:)
            b(j) = b(j) - f*b(i)
         end do
      end do
      do i = 5, 1, -1
         b(i) = (b(i) - sum(a(i,i+1:5)*b(i+1:5)))/a(i,i)
      end do
   end subroutine

   real(dp) function de2d(h, delta)
      ! int_0^1 dt int_0^{1-delta} du Ga(t,u)
      real(dp), intent(in) :: h, delta
      real(dp), allocatable :: xs(:), ws(:)
      real(dp) :: acc, ti, sc
      integer :: n, i, j
      call de_nodes(h, xs, ws, n)
      sc = 1 - delta
      acc = 0
      !$omp parallel do private(i,j,ti) reduction(+:acc)
      do i = 1, n
         ti = 0
         do j = 1, n
            ti = ti + ws(j)*gafun(xs(i), sc*xs(j))
         end do
         acc = acc + ws(i)*ti
      end do
      de2d = acc*sc
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
