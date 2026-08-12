! Diagram IIe: vacuum-polarization insertion, numeric check.
!
! mu_IIe = int_4^oo dt rho(t) K(t),  rho(t) = (1/3t)(1+2/t)sqrt(1-4/t),
! K(t)   = int_0^1 dz z(1-z)^2/((1-z)^2 + z t)          (m = 1)
!
! in units (alpha/pi)^2.  Substitution t = (1+y)^2/y maps t in (4,oo)
! to y in (0,1); the integrand is evaluated from t directly so this is
! an independent check of the SymPy transcription.
!
! Target: 119/36 - pi^2/3 = 0.0156874...
program g2_iie
   implicit none
   integer, parameter :: dp = kind(1.d0)
   integer, parameter :: n = 128
   real(dp) :: xg(n), wg(n)
   real(dp) :: s, y, z, t, dtdy, rho, kz, pi
   integer :: i, j

   pi = 4*atan(1._dp)
   call gauss_legendre_01(n, xg, wg)

   s = 0
   do i = 1, n
      y = xg(i)
      t = (1 + y)**2/y
      dtdy = (1 - y**2)/y**2          ! |dt/dy|
      rho = (1/(3*t))*(1 + 2/t)*sqrt(1 - 4/t)
      do j = 1, n
         z = xg(j)
         kz = z*(1 - z)**2/((1 - z)**2 + z*t)
         s = s + wg(i)*wg(j)*rho*dtdy*kz
      end do
   end do

   print "(a, f20.15)", "mu_IIe (numeric)      = ", s
   print "(a, f20.15)", "119/36 - pi^2/3       = ", 119._dp/36 - pi**2/3
   print "(a, es10.2)", "difference            = ", abs(s - (119._dp/36 - pi**2/3))

contains

   subroutine gauss_legendre_01(np, x, w)
      ! Gauss-Legendre nodes/weights on (0,1) by Newton iteration.
      integer, intent(in) :: np
      real(dp), intent(out) :: x(np), w(np)
      real(dp) :: xi, p0, p1, p2, dp1
      integer :: k, iter, l
      do k = 1, np
         xi = cos(pi*(k - 0.25_dp)/(np + 0.5_dp))
         do iter = 1, 100
            p0 = 1
            p1 = xi
            do l = 2, np
               p2 = ((2*l - 1)*xi*p1 - (l - 1)*p0)/l
               p0 = p1
               p1 = p2
            end do
            dp1 = np*(xi*p1 - p0)/(xi**2 - 1)
            if (abs(p1/dp1) < 1e-15_dp) exit
            xi = xi - p1/dp1
         end do
         x(k) = (1 - xi)/2
         w(k) = 1/((1 - xi**2)*dp1**2)
      end do
   end subroutine

end program
