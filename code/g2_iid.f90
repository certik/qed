! Diagram IId: numeric integration of the parametric integrands produced
! by g2_iid.py (g2_iid_frat.inc, g2_iid_flog.inc).
!
!   mu_IId(lam) = int f_rat dy dz du            (y+z<1; u in (0,1))
!               + int f_log dy dz dt du dxi     (y+z+t<1; u, xi in (0,1))
!
! Target (Petermann 1957, eq. (4)):
!   mu_IId = 11/24 - pi^2/18 + (1/2) log(lam^2)  =  -0.0899780... + log(lam)
!
! The program evaluates mu_IId for several photon masses lam and prints
! mu_IId - log(lam), which must approach -0.0899780 as lam -> 0.
module iid_integrands
   implicit none
   integer, parameter :: dp = kind(1.d0)
contains

   real(dp) function frat(y, z, u, lam)
      real(dp), intent(in) :: y, z, u, lam
      real(dp) :: f_rat
      include "g2_iid_frat.inc"
      frat = f_rat
   end function

   real(dp) function flog(y, z, t, u, xi, lam)
      real(dp), intent(in) :: y, z, t, u, xi, lam
      real(dp) :: f_log, C
      include "g2_iid_flog.inc"
      flog = f_log
   end function

end module

program g2_iid
   use iid_integrands
   implicit none
   integer, parameter :: nr = 240, nl = 32
   real(dp) :: xr(nr), wr(nr), xl(nl), wl(nl)
   real(dp) :: lam, s_rat, s_log, total, pi, target, c0
   real(dp) :: y, z, t, u, xi, jac, acc, vals(6)
   real(dp), parameter :: lams(6) = [0.1_dp, 0.03_dp, 0.01_dp, 0.003_dp, &
                                     0.001_dp, 0.0003_dp]
   integer :: i, j, a, b, c, il

   pi = 4*atan(1._dp)
   target = 11._dp/24 - pi**2/18
   call gauss_legendre_01(nr, xr, wr)
   call gauss_legendre_01(nl, xl, wl)
   ! smoothstep map u = 3s^2 - 2s^3 clusters points at both endpoints,
   ! resolving the IR structure at parameter-space scale lam
   do i = 1, nr
      wr(i) = wr(i)*6*xr(i)*(1 - xr(i))
      xr(i) = xr(i)**2*(3 - 2*xr(i))
   end do
   do i = 1, nl
      wl(i) = wl(i)*6*xl(i)*(1 - xl(i))
      xl(i) = xl(i)**2*(3 - 2*xl(i))
   end do

   print "(a)", "   lam     mu_IId - log(lam)    target = -0.0899780..."
   do il = 1, size(lams)
      lam = lams(il)

      ! rational piece: z in (0,1), y = (1-z) y', u in (0,1)
      s_rat = 0
      !$omp parallel do private(i,j,a,z,y,u,jac,acc) reduction(+:s_rat)
      do i = 1, nr
         z = xr(i)
         acc = 0
         do j = 1, nr
            y = (1 - z)*xr(j)
            jac = (1 - z)
            do a = 1, nr
               u = xr(a)
               acc = acc + wr(j)*wr(a)*jac*frat(y, z, u, lam)
            end do
         end do
         s_rat = s_rat + wr(i)*acc
      end do

      ! log piece: z, t = (1-z) t', y = (1-z)(1-t') y', u, xi
      s_log = 0
      !$omp parallel do private(i,j,a,b,c,z,t,y,u,xi,jac,acc) reduction(+:s_log)
      do i = 1, nl
         z = xl(i)
         acc = 0
         do j = 1, nl
            t = (1 - z)*xl(j)
            do a = 1, nl
               y = (1 - z)*(1 - xl(j))*xl(a)
               jac = (1 - z)**2*(1 - xl(j))
               do b = 1, nl
                  u = xl(b)
                  do c = 1, nl
                     xi = xl(c)
                     acc = acc + wl(j)*wl(a)*wl(b)*wl(c)*jac &
                           *flog(y, z, t, u, xi, lam)
                  end do
               end do
            end do
         end do
         s_log = s_log + wl(i)*acc
      end do

      total = s_rat + s_log
      vals(il) = total - log(lam)
      print "(f8.4, f20.12, a, f12.7, a, f12.7, a)", lam, vals(il), &
         "   (rat =", s_rat, ", log =", s_log, ")"
   end do
   ! extrapolate lam -> 0 with c0 + c1 lam + c2 lam log(lam) on the last 3
   call fit3(lams(4:6), vals(4:6), c0)
   print "(a, f20.12)", "extrapolated lam->0: ", c0
   print "(a, f20.12)", "target 11/24-pi^2/18:", target

contains

   subroutine fit3(ls, vs, c0)
      ! solve vs = c0 + c1 l + c2 l log(l) for the three points
      real(dp), intent(in) :: ls(3), vs(3)
      real(dp), intent(out) :: c0
      real(dp) :: m1(3, 3), r(3), det
      integer :: ii
      do ii = 1, 3
         m1(ii, :) = [1._dp, ls(ii), ls(ii)*log(ls(ii))]
         r(ii) = vs(ii)
      end do
      det = m1(1,1)*(m1(2,2)*m1(3,3) - m1(2,3)*m1(3,2)) &
          - m1(1,2)*(m1(2,1)*m1(3,3) - m1(2,3)*m1(3,1)) &
          + m1(1,3)*(m1(2,1)*m1(3,2) - m1(2,2)*m1(3,1))
      c0 = (r(1)*(m1(2,2)*m1(3,3) - m1(2,3)*m1(3,2)) &
          - m1(1,2)*(r(2)*m1(3,3) - m1(2,3)*r(3)) &
          + m1(1,3)*(r(2)*m1(3,2) - m1(2,2)*r(3)))/det
   end subroutine

   subroutine gauss_legendre_01(np, xq, wq)
      integer, intent(in) :: np
      real(dp), intent(out) :: xq(np), wq(np)
      real(dp) :: t0, p0, p1, p2, dp1
      integer :: kk, iter, l
      do kk = 1, np
         t0 = cos(pi*(kk - 0.25_dp)/(np + 0.5_dp))
         do iter = 1, 100
            p0 = 1
            p1 = t0
            do l = 2, np
               p2 = ((2*l - 1)*t0*p1 - (l - 1)*p0)/l
               p0 = p1
               p1 = p2
            end do
            dp1 = np*(t0*p1 - p0)/(t0**2 - 1)
            if (abs(p1/dp1) < 1e-15_dp) exit
            t0 = t0 - p1/dp1
         end do
         xq(kk) = (1 - t0)/2
         wq(kk) = 1/((1 - t0**2)*dp1**2)
      end do
   end subroutine

end program
