! Diagram I (crossed ladder) at lam = 0:
!    mu_I = int fI du dv dr dy dt,  u+v+r < 1,  y+t < 1
! evaluated by a randomly shifted rank-1 lattice rule (Korobov generating
! vector) over the 5-dim unit cube, in MIXED PRECISION.
!
! Each unit coordinate is passed through the order-3 smoothstep
!    s(x) = x^3 (10 - 15 x + 6 x^2),   s'(x) = 30 x^2 (1-x)^2,
! whose first two derivatives vanish at both ends: this concentrates
! points away from the edges and makes the periodic extension C^2, so the
! shifted lattice rule converges much faster than plain Monte Carlo.
! Independent random shifts give an honest error bar (each shift is an
! unbiased estimate regardless of the quality of the generating vector).
!
! fI itself is finite everywhere on the domain, but *evaluating* the
! generated expression in double precision loses all significance near
! the u-integration edges (u -> 0 and u+v+r -> 1), returning NaN within
! ~1e-5 of them; e.g. at (u, v, r, y, t) = (1e-8, .17, .82, .27, .73)
! double gives NaN while quad gives 0.101531754262515535.  Points closer
! than EDGE to either edge are therefore re-evaluated in real(16)
! (~3% of points, so the cost is dominated by the double-precision bulk).
! This needs gfortran: flang has no real(16) on macOS arm64.
!
! Product Gauss-Legendre failed here (boundary behaviour => the n-ladder
! is not monotone) and 5-dim tanh-sinh is far too expensive.
!
! target: 1/6 + 13/36 pi^2 + 5/4 zeta3 - 5/6 pi^2 log2 = -0.467645446094
module i_qmc_mod
   implicit none
   integer, parameter :: dp = kind(1.d0), qp = selected_real_kind(30)
   integer(8) :: nquad = 0, nbad = 0
   real(kind(1.d0)) :: wdrop = 0, wtot = 0   ! dropped vs total measure
   interface logc
      module procedure logc_r, logc_c, logc_rq, logc_cq
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
   complex(qp) function logc_rq(x)
      real(qp), intent(in) :: x
      logc_rq = log(cmplx(x, 0.0_qp, qp))
   end function
   complex(qp) function logc_cq(zz)
      complex(qp), intent(in) :: zz
      logc_cq = log(zz)
   end function

   real(dp) function ffI(u, v, r, y, t, lam)
      real(dp), intent(in) :: u, v, r, y, t, lam
      complex(dp) :: fIv
      complex(dp), parameter :: CI = (0.0_dp, 1.0_dp)
      ! g2_i_f_split.inc = g2_i_f.inc split into per-term statements by
      ! code/g2_i_split.py (compilers choke on the single 850 KB
      ! expression: >20 min and >30 GB for flang at any -O level)
      include "g2_i_f_split_decl.inc"
      include "g2_i_f_split.inc"
      ffI = real(fIv, dp)
   end function

   real(qp) function ffI_qp(u, v, r, y, t, lam)
      real(qp), intent(in) :: u, v, r, y, t, lam
      complex(qp) :: fIv
      complex(qp), parameter :: CI = (0.0_qp, 1.0_qp)
      ! same expression with quad literals (g2_i_split.py --quad)
      include "g2_i_f_split_qp_decl.inc"
      include "g2_i_f_split_qp.inc"
      ffI_qp = real(fIv, qp)
   end function

   real(dp) function feval(u, v, r, y, t, w)
      ! double precision in the bulk, quad where it breaks down.
      ! Calibration of double against quad on random points (worst
      ! relative error, code/g2_i_cal.f90):
      !    1-v-r = 1e-1: 1.8e-08     1e-3: 4.6e-01
      !            1e-2: 3.6e-04    <1e-4: 1.0e+00 (no digits left)
      ! so the fallback triggers on 1-v-r, not on u: the outer-loop
      ! b_hat = (r+v)(1-r-v) vanishing is what drives the cancellation.
      ! Below FLOOR even quad would lose its digits; those points carry
      ! a relative measure ~1e-6 under the smoothstep and a bounded
      ! integrand, and are dropped (counted in nbad).
      real(dp), intent(in) :: u, v, r, y, t, w
      real(dp), parameter :: EDGE = 1e-2_dp, FLOOR = 1e-6_dp
      real(dp) :: gv
      !$omp atomic
      wtot = wtot + w
      if (1 - v - r < FLOOR) then
         !$omp atomic
         nbad = nbad + 1
         !$omp atomic
         wdrop = wdrop + w
         feval = 0
         return
      else if (1 - v - r < EDGE) then
         !$omp atomic
         nquad = nquad + 1
         gv = real(ffI_qp(real(u, qp), real(v, qp), real(r, qp), &
                          real(y, qp), real(t, qp), 0.0_qp), dp)
      else
         gv = ffI(u, v, r, y, t, 0.0_dp)
      end if
      if (gv /= gv) then
         !$omp atomic
         nbad = nbad + 1
         !$omp atomic
         wdrop = wdrop + w
         gv = 0
      end if
      feval = gv
   end function
end module

program g2_i_qmc
   use i_qmc_mod
   implicit none
   integer(8), parameter :: agen = 1076671_8        ! Korobov parameter
   integer(8) :: nn = 16777216_8                    ! 2^24 lattice points
   integer :: nsh = 8                               ! random shifts
   real(dp) :: zv(5), sh(5), mean, sd, pi, targ
   real(dp), allocatable :: vals(:)
   real(dp), parameter :: zeta3 = 1.2020569031595942854_dp
   integer :: is, j, nargs
   integer(8) :: az
   integer, allocatable :: seed(:)
   character(len=8) :: argbuf
   ! optional arguments: log2(N) and the number of shifts, for quick runs
   nargs = command_argument_count()
   if (nargs >= 1) then
      call get_command_argument(1, argbuf)
      read (argbuf, *) j
      nn = 2_8**j
   end if
   if (nargs >= 2) then
      call get_command_argument(2, argbuf)
      read (argbuf, *) nsh
   end if
   allocate(vals(nsh))
   pi = 4*atan(1.0_dp)
   targ = 1.0_dp/6 + 13*pi**2/36 + 5*zeta3/4 - 5*pi**2*log(2.0_dp)/6
   call random_seed(size=j)
   allocate(seed(j))
   seed = 20260813
   call random_seed(put=seed)
   az = 1
   do j = 1, 5
      zv(j) = real(az, dp)
      az = mod(az*agen, nn)
   end do
   print "(a,i10,a,i3,a)", "lattice rule: N =", nn, ",", nsh, " shifts"
   do is = 1, nsh
      call random_number(sh)
      vals(is) = lattice_sum(zv, sh, nn)
      print "(a,i3,f18.10)", "  shift", is, vals(is)
      flush (6)
   end do
   mean = sum(vals)/nsh
   sd = sqrt(sum((vals - mean)**2)/(nsh*(nsh - 1.0_dp)))
   print "(a)", ""
   print "(a,f18.10,a,es9.2)", "mu_I   = ", mean, "  +/- ", sd
   print "(a,f18.10)", "target = ", targ
   print "(a,f18.10)", "diff   = ", mean - targ
   print "(a,f8.4,a)", "quad-precision evaluations: ", &
      100.0_dp*nquad/(nn*real(nsh, dp)), " %"
   print "(a,i12,a,es9.2)", "points dropped:             ", nbad, &
      "   fraction of the measure: ", wdrop/wtot
contains
   real(dp) function lattice_sum(zg, shift, np)
      real(dp), intent(in) :: zg(5), shift(5)
      integer(8), intent(in) :: np
      real(dp) :: acc, x(5), w, u, v, r, y, t, xr
      integer(8) :: i
      integer :: k
      acc = 0
      ! dynamic: the lattice's first coordinate is monotonic in i, so
      ! static chunks have wildly different quad-fallback loads
      !$omp parallel do private(i,k,x,w,u,v,r,y,t,xr) reduction(+:acc) &
      !$omp schedule(dynamic, 4096)
      do i = 0, np - 1
         w = 1
         do k = 1, 5
            xr = modulo(real(i, dp)*zg(k)/real(np, dp) + shift(k), 1.0_dp)
            x(k) = xr**3*(10 - 15*xr + 6*xr**2)      ! smoothstep
            w = w*30*xr**2*(1 - xr)**2               ! its derivative
         end do
         if (w == 0) cycle
         v = x(1)
         r = (1 - v)*x(2)
         u = (1 - v - r)*x(3)
         t = x(4)
         y = (1 - t)*x(5)
         w = w*(1 - v)*(1 - v - r)*(1 - t)
         acc = acc + w*feval(u, v, r, y, t, w)
      end do
      lattice_sum = acc/real(np, dp)
   end function
end program
