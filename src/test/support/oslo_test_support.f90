! Shared by the OSLO tests: the test problems, their reference solutions, the solver lists and the checks.
!
! A test calls check() for every condition it verifies and test_summary() at the end, which stops with a nonzero exit
! status when a check failed, so that ctest reports the test as failed.

module oslo_test_support
  use, intrinsic :: iso_fortran_env, only : R8 => real64, int64
  implicit none
  private
  public :: R8
  public :: check, test_summary
  public :: stiff_solvers, hairer_jacobian_solvers, solver_succeeded, error_ratio
  public :: robertson, robertson_jac, NROB, ROB_TIMES, ROB_REF
  public :: prothero_robinson, PR_LAMBDA, oscillator

  integer, save :: nchecks = 0, nfailed = 0

  ! Robertson's autocatalytic reaction (H. H. Robertson, The solution of a set of reaction rate equations, in J. Walsh
  ! (Ed.), Numerical Analysis: An Introduction, pp. 178-182, Academic Press, London, 1966):
  !    y1' = -0.04 y1 + 1e4 y2 y3
  !    y2' =  0.04 y1 - 1e4 y2 y3 - 3e7 y2^2
  !    y3' =  3e7 y2^2
  ! with y(0) = (1, 0, 0). The problem is stiff, and y1 + y2 + y3 = 1 for every t.
  ! The reference values at ROB_TIMES come from SciPy's Radau, BDF and LSODA with rtol = 1e-13, which agree to 3e-11
  ! (relative) at every time; they are given to 13 digits.
  integer, parameter :: NROB = 3
  real(R8), parameter :: ROB_TIMES(6) = [0.4_R8, 4.e1_R8, 4.e3_R8, 4.e5_R8, 4.e7_R8, 4.e10_R8]
  real(R8), parameter :: ROB_REF(NROB,6) = reshape([ &
    9.851721138610e-01_R8, 3.386395378975e-05_R8, 1.479402218522e-02_R8, &
    7.158270687194e-01_R8, 9.185534764558e-06_R8, 2.841637457458e-01_R8, &
    1.832022577767e-01_R8, 8.942371252776e-07_R8, 8.167968479862e-01_R8, &
    4.938274520980e-03_R8, 1.984994087954e-08_R8, 9.950617056291e-01_R8, &
    5.203071844121e-05_R8, 2.081335731893e-10_R8, 9.999479690734e-01_R8, &
    5.208345176798e-08_R8, 2.083338177925e-13_R8, 9.999999479163e-01_R8], [NROB,6])

  ! Prothero-Robinson: y' = lambda (y - sin t) + cos t, y(0) = 0, whose solution is sin t for any lambda; stiff for a
  ! large negative lambda
  real(R8), parameter :: PR_LAMBDA = -1.e4_R8

contains

  !> Records one check: a failed one is printed with its label
  subroutine check(condition, label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    nchecks = nchecks + 1
    if (.not. condition) then
      nfailed = nfailed + 1
      write(*,'(a)') 'FAILED: '//trim(label)
    end if
  end subroutine check

  !> Prints the count of the checks and stops with a nonzero exit status if one failed
  subroutine test_summary()
    write(*,'(/,i0,a,i0,a)') nchecks - nfailed, ' of ', nchecks, ' checks passed'
    if (nfailed > 0) error stop 1
  end subroutine test_summary

  !> The solvers for stiff problems that this build provides: all but the explicit H-dopri5 and, in the builds that
  !> disable them (ifx with OpenMP, see CMakeLists.txt), FATODE's Rosenbrock solvers
  function stiff_solvers() result(list)
    character(len=10), allocatable :: list(:)
    list = [character(len=10) :: 'H-radau5', 'radau2a', 'lobatto3c', 'gauss', 'radau1a', 'H-rodas', &
                                 'H-sdirk4', 'sdirk2a', 'sdirk2b', 'sdirk3a', 'sdirk4a', 'sdirk4b']
#   if !defined(FATODE_ROS_DISABLED)
    list = [character(len=10) :: list, 'ros2', 'ros3', 'ros4', 'rodas3', 'rodas4']
#   endif
#   if defined(SUNDIALS)
    list = [character(len=10) :: list, 'cvode']
#   endif
  end function stiff_solvers

  !> The solvers that use the analytical Jacobian passed to run_odesolver (IJAC = 1)
  function hairer_jacobian_solvers() result(list)
    character(len=10), allocatable :: list(:)
    list = [character(len=10) :: 'H-radau5', 'H-rodas', 'H-sdirk4']
  end function hairer_jacobian_solvers

  !> Whether the status that run_odesolver returned means success. Each family keeps its own convention: Hairer's
  !> IDID and FATODE's IERR are positive on success, CVODE's flag is non-negative
  logical function solver_succeeded(solver, err)
    character(len=*), intent(in) :: solver
    integer, intent(in) :: err
    select case (solver)
    case ('cvode')
      solver_succeeded = err >= 0
    case default
      solver_succeeded = err > 0
    end select
  end function solver_succeeded

  !> The largest error of y against yref in units of the tolerance, |y - yref| / (rt |yref| + at); a NaN or an
  !> infinity in y gives a huge value, so that it fails any bound
  real(R8) function error_ratio(y, yref, rt, at)
    real(R8), intent(in) :: y(:), yref(:), rt(:), at(:)
    error_ratio = maxval(abs(y - yref) / (rt*abs(yref) + at))
    if (any(not_finite(y))) error_ratio = huge(1._R8)
  end function error_ratio

  !> Whether x is a NaN or an infinity, from its exponent bits: the release flags (-ffast-math, -fp-model fast) let the
  !> compiler assume that no value is a NaN and drop a test such as x /= x
  elemental logical function not_finite(x)
    real(R8), intent(in) :: x
    not_finite = iand(shiftr(transfer(x, 0_int64), 52), 2047_int64) == 2047_int64
  end function not_finite

  subroutine robertson(neq, t, y, ydot)
    integer, intent(in) :: neq
    real(R8), intent(in) :: t
    real(R8), intent(in) :: y(neq)
    real(R8), intent(out) :: ydot(neq)
    ydot(1) = -0.04_R8*y(1) + 1.e4_R8*y(2)*y(3)
    ydot(2) = 0.04_R8*y(1) - 1.e4_R8*y(2)*y(3) - 3.e7_R8*y(2)**2
    ydot(3) = 3.e7_R8*y(2)**2
  end subroutine robertson

  subroutine robertson_jac(n, x, y, dfy, ldfy, rpar, ipar)
    integer, intent(in) :: n, ldfy
    real(R8), intent(in) :: x, y(n)
    real(R8), intent(out) :: dfy(ldfy,n)
    real(R8), intent(in) :: rpar(*)
    integer, intent(in) :: ipar(*)
    dfy(1,:) = [-0.04_R8,  1.e4_R8*y(3),                    1.e4_R8*y(2)]
    dfy(2,:) = [ 0.04_R8, -1.e4_R8*y(3) - 6.e7_R8*y(2), -1.e4_R8*y(2)]
    dfy(3,:) = [ 0._R8,    6.e7_R8*y(2),                    0._R8]
  end subroutine robertson_jac

  subroutine prothero_robinson(neq, t, y, ydot)
    integer, intent(in) :: neq
    real(R8), intent(in) :: t
    real(R8), intent(in) :: y(neq)
    real(R8), intent(out) :: ydot(neq)
    ydot(1) = PR_LAMBDA*(y(1) - sin(t)) + cos(t)
  end subroutine prothero_robinson

  !> The harmonic oscillator y1' = y2, y2' = -y1, whose solution from y(0) = (1, 0) is (cos t, -sin t): not stiff
  subroutine oscillator(neq, t, y, ydot)
    integer, intent(in) :: neq
    real(R8), intent(in) :: t
    real(R8), intent(in) :: y(neq)
    real(R8), intent(out) :: ydot(neq)
    ydot(1) = y(2)
    ydot(2) = -y(1)
  end subroutine oscillator

end module oslo_test_support
