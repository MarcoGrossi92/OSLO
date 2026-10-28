! Accuracy of every solver on two problems with a known solution (see oslo_test_support), checked at four times:
!   - Prothero-Robinson with lambda = -1e4 (stiff), solution sin t, for every stiff solver;
!   - the harmonic oscillator (not stiff), solution (cos t, -sin t), for every solver, the explicit H-dopri5 included,
!     but radau1a: FATODE's error estimators do not control Radau IA on it. With the SDIRK estimator, which OSLO uses,
!     radau1a accepts steps of nearly a period and ends with an error of order 1; with the classic one it is accurate
!     there but fails the same way on Prothero-Robinson.
! Every call must return a successful status and stay within a multiple of the tolerance (rt |y| + at) of the exact
! solution. The bounds are an order of magnitude above what the solvers reach with gfortran; the second-order methods
! and CVODE's BDF, whose global error on the oscillator grows to 20-50 times the tolerance, have their own.

program exact_solutions_test
  use oslo
  use oslo_test_support
  implicit none
  real(R8), parameter :: RT = 1.e-8_R8, AT = 1.e-10_R8
  real(R8), parameter :: MAX_RATIO = 50._R8, MAX_RATIO_LOOSE = 250._R8
  real(R8), parameter :: TIMES(4) = [2.5_R8, 5._R8, 7.5_R8, 10._R8]
  character(len=10), allocatable :: solvers(:)
  integer :: i

  write(*,'(a10,a20,a14)') 'solver', 'problem', 'max err/tol'

  solvers = stiff_solvers()
  do i = 1, size(solvers)
    call run_exact(trim(solvers(i)), 'prothero-robinson')
  end do

  solvers = [character(len=10) :: 'H-dopri5', stiff_solvers()]
  do i = 1, size(solvers)
    if (solvers(i) == 'radau1a') cycle
    call run_exact(trim(solvers(i)), 'oscillator')
  end do

  call test_summary()

contains

  subroutine run_exact(solver, problem)
    character(len=*), intent(in) :: solver, problem
    real(R8), allocatable :: y(:), yexact(:)
    real(R8) :: t, tout, ratio, max_ratio_seen, ratio_bound
    integer :: k, n, err
    logical :: ok
    character(len=64) :: label

    ratio_bound = MAX_RATIO
    if (any(solver == [character(len=8) :: 'ros2', 'sdirk2a', 'sdirk2b', 'cvode'])) ratio_bound = MAX_RATIO_LOOSE

    if (problem == 'oscillator') then
      n = 2
      y = [1._R8, 0._R8]
    else
      n = 1
      y = [0._R8]
    end if
    call setup_odesolver(N=n, solver=solver, RT=spread(RT, 1, n), AT=spread(AT, 1, n))
    t = 0._R8
    max_ratio_seen = 0._R8
    ok = .true.
    do k = 1, size(TIMES)
      tout = TIMES(k)
      if (problem == 'oscillator') then
        call run_odesolver(n, t, tout, y, oscillator, no_jacobian, 0, err)
        yexact = [cos(TIMES(k)), -sin(TIMES(k))]
      else
        call run_odesolver(n, t, tout, y, prothero_robinson, no_jacobian, 0, err)
        yexact = [sin(TIMES(k))]
      end if
      write(label,'(a,1x,a,a,f4.1)') solver, problem, ' at t = ', TIMES(k)
      ok = solver_succeeded(solver, err)
      write(label(len_trim(label)+1:),'(a,i0)') ', status ', err
      call check(ok, trim(label)//': unsuccessful')
      if (.not. ok) exit
      ratio = error_ratio(y, yexact, spread(RT, 1, n), spread(AT, 1, n))
      call check(ratio <= ratio_bound, trim(label)//': error above the bound')
      max_ratio_seen = max(max_ratio_seen, ratio)
      t = TIMES(k)
    end do

    if (ok) then
      write(*,'(a10,a20,es14.3)') solver, problem, max_ratio_seen
    else
      write(*,'(a10,a20,a)') solver, problem, '        failed'
    end if
  end subroutine run_exact

end program exact_solutions_test
