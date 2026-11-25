! The perturbation SDIRK4 gives a variable for its finite-difference Jacobian: sqrt(u max(1e-5,|y|)) with u the rounding
! unit of double precision, 1.1e-16. With y = 1 it is 1.05e-8. A smaller u amplifies the rounding of the right-hand side.

module jacobian_probe
  use oslo_test_support, only : R8
  implicit none
  real(R8) :: smallest
contains
  subroutine probe_rhs(neq, t, y, ydot)
    integer, intent(in) :: neq
    real(R8), intent(in) :: t, y(neq)
    real(R8), intent(out) :: ydot(neq)
    ! y(2) has ydot = 0: it stays at 1, and only the Jacobian moves it
    if (y(2) /= 1._R8) smallest = min(smallest, abs(y(2) - 1._R8))
    ydot(1) = -1.e3_R8*(y(1) - cos(t))
    ydot(2) = 0._R8
  end subroutine probe_rhs
end module jacobian_probe

program sdirk4_jacobian_step
  use oslo
  use oslo_test_support
  use jacobian_probe
  implicit none
  real(R8) :: y(2), t1, t2, rt(2), at(2)
  integer :: err, io(3)
  rt = 1.e-8_R8; at = 1.e-10_R8; io = 0
  smallest = huge(1._R8)
  y = [0._R8, 1._R8]; t1 = 0._R8; t2 = 1._R8
  call setup_odesolver(N=2, solver='H-sdirk4', RT=rt, AT=at, iopt=io)
  call run_odesolver(2, t1, t2, y, probe_rhs, no_jacobian, 0, err)
  print '(a,es10.3)', ' smallest perturbation of a variable equal to 1: ', smallest
  if (.not. (err > 0)) error stop 'H-sdirk4 failed'
  if (smallest < 0.5e-8_R8 .or. smallest > 2.e-8_R8) error stop 'the Jacobian perturbation is not that of double precision'
end program sdirk4_jacobian_step
