! Accuracy of every stiff solver on Robertson's problem (see oslo_test_support): the solution is checked at six times
! from 0.4 to 4e10 against reference values computed independently of OSLO. Every solver must
!   - return a successful status from each call,
!   - stay within a multiple of the tolerance (rt |y| + at) of the reference, for each component and time,
!   - keep y1 + y2 + y3 = 1, which the exact solution and every method in exact arithmetic do.
! The bounds are an order of magnitude above what the solvers reach with gfortran. They are looser for the Rosenbrock
! methods with a finite-difference Jacobian: their stages use the Jacobian itself, not only to iterate Newton, so its
! error goes into the solution (with the analytical Jacobian they are as accurate as the others).
! The solvers that accept an analytical Jacobian (IJAC = 1) are run with it as well.

program robertson_test
  use oslo
  use oslo_test_support
  implicit none
  real(R8), parameter :: RT(NROB) = 1.e-6_R8
  real(R8), parameter :: AT(NROB) = [1.e-12_R8, 1.e-18_R8, 1.e-12_R8]
  real(R8), parameter :: MAX_RATIO = 50._R8, MAX_RATIO_ROS_FD = 2000._R8
  real(R8), parameter :: INVARIANT_TOL = 1.e-10_R8, INVARIANT_TOL_ROS_FD = 1.e-7_R8
  character(len=10), allocatable :: solvers(:)
  integer :: i

  write(*,'(a10,a6,a14,a14)') 'solver', 'ijac', 'max err/tol', '|sum(y) - 1|'

  solvers = stiff_solvers()
  do i = 1, size(solvers)
    call run_robertson(trim(solvers(i)), 0)
  end do

  solvers = hairer_jacobian_solvers()
  do i = 1, size(solvers)
    call run_robertson(trim(solvers(i)), 1)
  end do

  call test_summary()

contains

  subroutine run_robertson(solver, ijac)
    character(len=*), intent(in) :: solver
    integer, intent(in) :: ijac
    real(R8) :: y(NROB), t, tout, ratio, max_ratio_seen, drift, max_drift, ratio_bound, drift_bound
    integer :: k, err
    logical :: ok
    character(len=64) :: label

    ratio_bound = MAX_RATIO
    drift_bound = INVARIANT_TOL
    if (ijac == 0 .and. is_rosenbrock(solver)) then
      ratio_bound = MAX_RATIO_ROS_FD
      drift_bound = INVARIANT_TOL_ROS_FD
    end if

    call setup_odesolver(N=NROB, solver=solver, RT=RT, AT=AT)
    y = [1._R8, 0._R8, 0._R8]
    t = 0._R8
    max_ratio_seen = 0._R8
    max_drift = 0._R8
    ok = .true.
    do k = 1, size(ROB_TIMES)
      tout = ROB_TIMES(k)
      if (ijac == 1) then
        call run_odesolver(NROB, t, tout, y, robertson, robertson_jac, 1, err)
      else
        call run_odesolver(NROB, t, tout, y, robertson, no_jacobian, 0, err)
      end if
      write(label,'(a,a,i0,a,es8.1)') solver, ' (ijac=', ijac, ') at t = ', ROB_TIMES(k)
      ok = solver_succeeded(solver, err)
      call check(ok, trim(label)//': status '//status_string(err))
      if (.not. ok) exit
      ratio = error_ratio(y, ROB_REF(:,k), RT, AT)
      drift = abs(sum(y) - 1._R8)
      call check(ratio <= ratio_bound, trim(label)//': error above the bound')
      call check(drift <= drift_bound, trim(label)//': y1 + y2 + y3 drifted from 1')
      max_ratio_seen = max(max_ratio_seen, ratio)
      max_drift = max(max_drift, drift)
      t = ROB_TIMES(k)
    end do

    if (ok) then
      write(*,'(a10,i6,2es14.3)') solver, ijac, max_ratio_seen, max_drift
    else
      write(*,'(a10,i6,a)') solver, ijac, '        failed'
    end if
  end subroutine run_robertson

  logical function is_rosenbrock(solver)
    character(len=*), intent(in) :: solver
    is_rosenbrock = any(solver == [character(len=8) :: 'H-rodas', 'ros2', 'ros3', 'ros4', 'rodas3', 'rodas4'])
  end function is_rosenbrock

  function status_string(err) result(s)
    integer, intent(in) :: err
    character(len=12) :: s
    write(s,'(i0)') err
  end function status_string

end program robertson_test
