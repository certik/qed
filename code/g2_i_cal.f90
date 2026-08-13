! Calibration for diagram I: where does double-precision evaluation of
! the generated fI lose its significance?  Compares ffI (double) against
! ffI_qp (quad) on random points at a ladder of distances from the
! 1-v-r = 0 face, and times both.  Used to set the EDGE threshold of the
! mixed-precision fallback in g2_i_qmc.f90.  Needs gfortran (real(16)).
module m
   implicit none
   integer, parameter :: dp = kind(1.d0), qp = selected_real_kind(30)
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
      include "g2_i_f_split_decl.inc"
      include "g2_i_f_split.inc"
      ffI = real(fIv, dp)
   end function
   real(qp) function ffI_qp(u, v, r, y, t, lam)
      real(qp), intent(in) :: u, v, r, y, t, lam
      complex(qp) :: fIv
      complex(qp), parameter :: CI = (0.0_qp, 1.0_qp)
      include "g2_i_f_split_qp_decl.inc"
      include "g2_i_f_split_qp.inc"
      ffI_qp = real(fIv, qp)
   end function
end module
program cal
   use m
   implicit none
   real(dp) :: u, v, r, y, t, gd, worst(9), t0, t1
   real(qp) :: gq
   real(dp) :: scales(9) = [1e-1_dp, 1e-2_dp, 1e-3_dp, 1e-4_dp, 1e-5_dp, &
                            1e-6_dp, 1e-7_dp, 1e-8_dp, 1e-9_dp]
   real(dp) :: rn(4)
   integer :: k, i, nrep
   integer, allocatable :: seed(:)
   call random_seed(size=k); allocate(seed(k)); seed = 7; call random_seed(put=seed)
   worst = 0
   print "(a)", " scale of (1-v-r)      :  worst relative error of double"
   do k = 1, 9
      do i = 1, 40
         call random_number(rn)
         v = 0.05_dp + 0.9_dp*rn(1)
         r = (1 - v) - scales(k)*(0.3_dp + 0.7_dp*rn(2))
         if (r <= 0) cycle
         u = (1 - v - r)*(0.2_dp + 0.6_dp*rn(3))
         t = 0.05_dp + 0.9_dp*rn(4)
         y = (1 - t)*0.4_dp
         gd = ffI(u, v, r, y, t, 0.0_dp)
         gq = ffI_qp(real(u, qp), real(v, qp), real(r, qp), real(y, qp), &
                     real(t, qp), 0.0_qp)
         if (gd /= gd) then
            worst(k) = 1
         else if (abs(gq) > 0) then
            worst(k) = max(worst(k), abs((real(gq, dp) - gd)/real(gq, dp)))
         end if
      end do
      print "(es10.1,es14.2)", scales(k), worst(k)
   end do
   nrep = 60
   call cpu_time(t0)
   do i = 1, nrep
      gq = ffI_qp(0.11_qp + i*1e-4_qp, 0.23_qp, 0.31_qp, 0.17_qp, 0.29_qp, 0.0_qp)
   end do
   call cpu_time(t1)
   print "(a,f10.1,a,es14.6)", "quad per eval: ", (t1 - t0)/nrep*1e6, " us", real(gq, dp)
   call cpu_time(t0)
   do i = 1, 20000
      gd = ffI(0.11_dp + i*1e-7_dp, 0.23_dp, 0.31_dp, 0.17_dp, 0.29_dp, 0.0_dp)
   end do
   call cpu_time(t1)
   print "(a,f10.1,a,es14.6)", "dble per eval: ", (t1 - t0)/20000*1e6, " us", gd
end program
