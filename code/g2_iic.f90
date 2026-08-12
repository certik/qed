! Diagram IIc (corner): numeric integration of the parametric integrands
! produced by g2_iic.py (g2_iic_fa.inc, g2_iic_fb.inc, g2_iic_fc.inc).
!
!   mu_IIc(lam) = int fa dy dz dt du dv       (y+z+t<1, u+v<1)
!               + int fb dy dz du dv          (y+z<1,   u+v<1)
!               + int fc dy dz dt du dv dxi   (y+z+t<1, u+v<1, xi in (0,1))
!
! Target (Petermann 1957, eq. (3)):
!   mu_IIc = -67/24 + pi^2/18 - zeta(3)/2 + (pi^2/3) log 2 - log(lam)
!          = -0.564021... - log(lam)
!
! The program prints mu_IIc + log(lam), which must -> -0.564021 as lam -> 0.
module iic_integrands
   implicit none
   integer, parameter :: dp = kind(1.d0)
contains

   real(dp) function ffa(y, z, t, u, v, lam)
      real(dp), intent(in) :: y, z, t, u, v, lam
      real(dp) :: fa
      include "g2_iic_fa.inc"
      ffa = fa
   end function

   real(dp) function ffb(y, z, u, v, lam)
      real(dp), intent(in) :: y, z, u, v, lam
      real(dp) :: fb
      include "g2_iic_fb.inc"
      ffb = fb
   end function

   real(dp) function ffc(y, z, t, u, v, xi, lam)
      real(dp), intent(in) :: y, z, t, u, v, xi, lam
      real(dp) :: fc
      include "g2_iic_fc.inc"
      ffc = fc
   end function

end module

program g2_iic
   use iic_integrands
   implicit none
   integer, parameter :: na = 28, nb = 96, nc = 18
   real(dp) :: xa(na), wa(na), xb(nb), wb(nb), xc(nc), wc(nc)
   real(dp) :: lam, sa, sb, sc, total, pi, target, c0
   real(dp) :: y, z, t, u, v, xi, jac, acc, fcv, vals(6)
   real(dp), parameter :: lams(6) = [0.1_dp, 0.03_dp, 0.01_dp, 0.003_dp, &
                                     0.001_dp, 0.0003_dp]
   real(dp), parameter :: zeta3 = 1.2020569031595942854_dp
   integer :: i, j, a, b, c, d, il

   pi = 4*atan(1._dp)
   target = -67._dp/24 + pi**2/18 - zeta3/2 + pi**2*log(2._dp)/3
   call gl01(na, xa, wa)
   call gl01(nb, xb, wb)
   call gl01(nc, xc, wc)
   call smoothstep(na, xa, wa)
   call smoothstep(nb, xb, wb)
   call smoothstep(nc, xc, wc)

   print "(a)", "   lam     mu_IIc + log(lam)    target = -0.5640214..."
   do il = 1, size(lams)
      lam = lams(il)

      ! fb: y+z<1 (y=(1-z)y'), u+v<1 (u=(1-v)u')
      sb = 0
      !$omp parallel do private(i,j,a,b,z,y,v,u,jac,acc) reduction(+:sb)
      do i = 1, nb
         z = xb(i)
         acc = 0
         do j = 1, nb
            y = (1 - z)*xb(j)
            do a = 1, nb
               v = xb(a)
               do b = 1, nb
                  u = (1 - v)*xb(b)
                  jac = (1 - z)*(1 - v)
                  acc = acc + wb(j)*wb(a)*wb(b)*jac*ffb(y, z, u, v, lam)
               end do
            end do
         end do
         sb = sb + wb(i)*acc
      end do

      ! fa: z, t=(1-z)t', y=(1-z)(1-t')y'; u+v<1
      sa = 0
      !$omp parallel do private(i,j,a,b,c,z,t,y,v,u,jac,acc) reduction(+:sa)
      do i = 1, na
         z = xa(i)
         acc = 0
         do j = 1, na
            t = (1 - z)*xa(j)
            do a = 1, na
               y = (1 - z)*(1 - xa(j))*xa(a)
               do b = 1, na
                  v = xa(b)
                  do c = 1, na
                     u = (1 - v)*xa(c)
                     jac = (1 - z)**2*(1 - xa(j))*(1 - v)
                     acc = acc + wa(j)*wa(a)*wa(b)*wa(c)*jac &
                           *ffa(y, z, t, u, v, lam)
                  end do
               end do
            end do
         end do
         sa = sa + wa(i)*acc
      end do

      ! fc: as fa plus xi
      sc = 0
      !$omp parallel do private(i,j,a,b,c,d,z,t,y,v,u,xi,jac,acc,fcv) reduction(+:sc)
      do i = 1, nc
         z = xc(i)
         acc = 0
         do j = 1, nc
            t = (1 - z)*xc(j)
            do a = 1, nc
               y = (1 - z)*(1 - xc(j))*xc(a)
               do b = 1, nc
                  v = xc(b)
                  do c = 1, nc
                     u = (1 - v)*xc(c)
                     do d = 1, nc
                        xi = xc(d)
                        jac = (1 - z)**2*(1 - xc(j))*(1 - v)
                        ! bhat = u(1-u) -> 0 gives numerical 0/0
                        ! (finite true limit); skip the sliver
                        if (u < 2e-4_dp .or. 1 - u < 2e-4_dp) then
                           fcv = 0
                        else
                           fcv = ffc(y, z, t, u, v, xi, lam)
                        end if
                        acc = acc + wc(j)*wc(a)*wc(b)*wc(c)*wc(d)*jac*fcv
                     end do
                  end do
               end do
            end do
         end do
         sc = sc + wc(i)*acc
      end do

      total = sa + sb + sc
      vals(il) = total + log(lam)
      print "(f8.4, f20.12, 3(a, f13.7))", lam, vals(il), &
         "   a=", sa, " b=", sb, " c=", sc
   end do
   call fit3(lams(4:6), vals(4:6), c0)
   print "(a, f20.12)", "extrapolated lam->0: ", c0
   print "(a, f20.12)", "target:              ", target

contains

   subroutine smoothstep(np, xq, wq)
      integer, intent(in) :: np
      real(dp), intent(inout) :: xq(np), wq(np)
      integer :: kk
      do kk = 1, np
         wq(kk) = wq(kk)*6*xq(kk)*(1 - xq(kk))
         xq(kk) = xq(kk)**2*(3 - 2*xq(kk))
      end do
   end subroutine

   subroutine fit3(ls, vs, c0)
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
