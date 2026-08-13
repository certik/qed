! Diagram I (crossed ladder): numeric lambda sweep of fI
! mu_I = int fI du dv dr dy dt  (u+v+r<1, y+t<1)
! target: 1/6 + 13/36 pi^2 + 5/4 zeta3 - 5/6 pi^2 log2 = -0.4676445...
module i_integrands
   implicit none
   integer, parameter :: dp = kind(1.d0)
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
      include "g2_i_f.inc"
      ffI = real(fIv, dp)
   end function
end module

program g2_i
   use i_integrands
   implicit none
   integer, parameter :: n = 24
   real(dp) :: xg(n), wg(n), lam, acc, pi, target
   real(dp) :: u, v, r, y, t, jac
   real(dp), parameter :: lams(4) = [0.1_dp, 0.03_dp, 0.01_dp, 0.003_dp]
   real(dp), parameter :: zeta3 = 1.2020569031595942854_dp
   real(dp) :: epsg
   integer :: nbad = 0
   integer :: i1, i2, i3, i4, i5, il
   pi = 4*atan(1.0_dp)
   target = 1.0_dp/6 + 13*pi**2/36 + 5*zeta3/4 - 5*pi**2*log(2.0_dp)/6
   call gl01(n, xg, wg)
   do i1 = 1, n
      wg(i1) = wg(i1)*6*xg(i1)*(1 - xg(i1))
      xg(i1) = xg(i1)**2*(3 - 2*xg(i1))
   end do
   print "(a)", "   lam     mu_I               target = -0.4676445..."
   do il = 1, size(lams)
      epsg = 1e-5_dp
      if (il == 1) epsg = 3e-5_dp   ! sensitivity check on first lam
      lam = lams(il)
      acc = 0
      !$omp parallel do private(i1,i2,i3,i4,i5,u,v,r,y,t,jac) reduction(+:acc)
      do i1 = 1, n
         v = xg(i1)
         do i2 = 1, n
            r = (1 - v)*xg(i2)
            do i3 = 1, n
               u = (1 - v - r)*xg(i3)
               do i4 = 1, n
                  t = xg(i4)
                  do i5 = 1, n
                     y = (1 - t)*xg(i5)
                     jac = (1 - v)*(1 - v - r)*(1 - t)
                     ! double-precision evaluation garbage near edges
                     if (min(u, v, r, y, t, 1 - u - v - r, 1 - y - t) &
                         > epsg) then
                        block
                           real(dp) :: gv
                           gv = ffI(u, v, r, y, t, lam)
                           if (gv /= gv) then
                              !$omp critical
                              nbad = nbad + 1
                              if (nbad <= 5) print "(a,5es11.2)", "bad:", &
                                 u, v, r, y, t
                              !$omp end critical
                              gv = 0
                           end if
                           acc = acc + wg(i1)*wg(i2)*wg(i3)*wg(i4)*wg(i5) &
                                 *jac*gv
                        end block
                     end if
                  end do
               end do
            end do
         end do
      end do
      print "(f8.4, f20.12)", lam, acc
   end do
   print "(a, f20.12)", "target: ", target
contains
   subroutine gl01(np, xq, wq)
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
