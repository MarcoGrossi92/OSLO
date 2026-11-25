! A right-hand side that has no value at a trial state (a NaN there) makes SDIRK4 reject the step and halve it: the
! solution is never a NaN and stays accurate. ydot = -K sqrt(y) + s has its NaN at y < 0, which the extrapolated
! starting values of a large step reach; its exact solution decays from 1 to the equilibrium (s/K)^2 and stays positive.
!   - any run: no NaN is returned, and success is only reported with the right result;
!   - a run whose equilibrium is at least 50 times the tolerance succeeds. Below that the error control may accept a
!     state on the wrong side of zero, where the problem has no right-hand side: such a run may end as a failure.

module trial_nan
  use oslo_test_support, only : R8
  implicit none
  real(R8) :: K, s
  integer :: nans
contains
  subroutine sqrt_rhs(neq, t, y, ydot)
    integer, intent(in) :: neq
    real(R8), intent(in) :: t, y(neq)
    real(R8), intent(out) :: ydot(neq)
    real(R8), volatile :: zero
    zero = 0._R8
    if (y(1) < 0._R8) then
      nans = nans + 1
      ydot(1) = zero/zero
    else
      ydot(1) = -K*sqrt(y(1)) + s
    end if
  end subroutine sqrt_rhs
end module trial_nan

program sdirk4_trial_nan
  use oslo
  use oslo_test_support
  use trial_nan
  implicit none
  real(R8), parameter :: KS(4) = [1.e1_R8, 1.e2_R8, 1.e3_R8, 1.e4_R8], SS(4) = [1.e-3_R8, 1.e-1_R8, 1._R8, 1.e1_R8]
  real(R8), parameter :: RTS(3) = [1.e-4_R8, 1.e-6_R8, 1.e-8_R8]
  real(R8) :: y(1), t1, t2, rt(1), at(1), yeq
  integer :: err, io(3), i, j, m, bad, total, required
  logical :: wrong, nan
  bad = 0; total = 0; required = 0; io = 0
  do i = 1, size(KS)
    do j = 1, size(SS)
      do m = 1, size(RTS)
        K = KS(i); s = SS(j); rt = RTS(m); at = RTS(m)
        nans = 0
        y = 1._R8; t1 = 0._R8; t2 = 10._R8
        call setup_odesolver(N=1, solver='H-sdirk4', RT=rt, AT=at, iopt=io)
        call run_odesolver(1, t1, t2, y, sqrt_rhs, no_jacobian, 0, err)
        yeq = (s/K)**2
        total = total + nans
        wrong = error_ratio(y, [yeq], rt, at) > 50._R8
        nan = error_ratio(y, [yeq], rt, at) == huge(1._R8)
        if (yeq >= 50._R8*at(1)) required = required + min(nans, 1)
        if (nan .or. (err > 0 .and. wrong) .or. (yeq >= 50._R8*at(1) .and. (err <= 0 .or. wrong))) then
          bad = bad + 1
          print '(a,3es9.1,a,i0,a,i0,a,es12.4,a,es12.4)', ' K s tol =', K, s, rt, '  err = ', err, '  NaN calls = ', nans, &
            '  y = ', y(1), '  equilibrium ', yeq
        end if
      end do
    end do
  end do
  print '(a,i0,a,i0)', ' NaN calls in all: ', total, '   wrong results: ', bad
  print '(a,i0)', ' runs with a NaN whose success is required: ', required
  if (required == 0) error stop 'no run that must succeed reaches a NaN: the test sees nothing'
  if (bad > 0) error stop 'wrong result'
end program sdirk4_trial_nan
