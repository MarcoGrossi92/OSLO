! The optional arguments of run_odesolver and setup_odesolver:
!   - solout (H-dopri5, H-sdirk4): called once at the start and after each accepted step, with consecutive step numbers,
!     steps that join up, the start and end times, and the solution of the step (the last one is what the solver returns);
!   - hmax: no step is longer;
!   - a solout that sets IRTRN < 0 stops the integration, which returns status 2 at the time it stopped, with the
!     solution there;
!   - the step limit IOPT(1): a solver that needs more steps returns a failure, which shows that the solvers do report
!     one (cvode does not take the limit from IOPT).

module solout_recorder
  use oslo_test_support, only : R8
  implicit none
  integer :: ncalls
  logical :: numbers_consecutive, steps_join
  real(R8) :: first_x, last_x, last_y(2), max_step, stop_at

contains

  subroutine reset_recorder(stop_time)
    real(R8), intent(in) :: stop_time
    ncalls = 0
    numbers_consecutive = .true.
    steps_join = .true.
    first_x = -huge(1._R8)
    last_x = -huge(1._R8)
    last_y = 0._R8
    max_step = 0._R8
    stop_at = stop_time
  end subroutine reset_recorder

  subroutine record_solout(NR, XOLD, X, Y, N, IRTRN)
    integer  :: NR, N
    real(R8) :: X, XOLD, Y(N)
    integer  :: IRTRN
    ncalls = ncalls + 1
    if (NR /= ncalls) numbers_consecutive = .false.
    if (ncalls == 1) then
      first_x = X
    else
      if (XOLD /= last_x .or. X <= XOLD) steps_join = .false.
      max_step = max(max_step, X - XOLD)
    end if
    last_x = X
    last_y = Y(1:2)
    if (X >= stop_at) IRTRN = -1
  end subroutine record_solout

end module solout_recorder


program options_test
  use oslo
  use oslo_test_support
  use solout_recorder
  implicit none
  real(R8), parameter :: RT = 1.e-8_R8, AT = 1.e-10_R8
  real(R8), parameter :: MAX_RATIO = 50._R8
  real(R8), parameter :: TEND = 10._R8, HMAX = 0.05_R8, STOP_TIME = 5._R8
  character(len=10), allocatable :: solvers(:)
  integer :: i

  solvers = [character(len=10) :: 'H-dopri5', 'H-sdirk4']
  do i = 1, size(solvers)
    call check_solout(trim(solvers(i)))
    call check_hmax(trim(solvers(i)))
    call check_interrupt(trim(solvers(i)))
  end do

  call check_step_limit('H-dopri5')
  solvers = stiff_solvers()
  do i = 1, size(solvers)
    if (solvers(i) == 'cvode') cycle
    call check_step_limit(trim(solvers(i)))
  end do

  call test_summary()

contains

  !> Integrates the oscillator from 0 to TEND, recording the calls of solout
  subroutine integrate_oscillator(solver, stop_time, y, t, err, hmax_arg)
    character(len=*), intent(in) :: solver
    real(R8), intent(in) :: stop_time
    real(R8), intent(out) :: y(2), t
    integer, intent(out) :: err
    real(R8), intent(in), optional :: hmax_arg
    real(R8) :: tout
    call setup_odesolver(N=2, solver=solver, RT=[RT, RT], AT=[AT, AT])
    call reset_recorder(stop_time)
    y = [1._R8, 0._R8]
    t = 0._R8
    tout = TEND
    if (present(hmax_arg)) then
      call run_odesolver(2, t, tout, y, oscillator, no_jacobian, 0, err, hmax=hmax_arg, solout=record_solout)
    else
      call run_odesolver(2, t, tout, y, oscillator, no_jacobian, 0, err, solout=record_solout)
    end if
  end subroutine integrate_oscillator

  subroutine check_solout(solver)
    character(len=*), intent(in) :: solver
    real(R8) :: y(2), t
    integer :: err
    call integrate_oscillator(solver, huge(1._R8), y, t, err)
    write(*,'(a10,a,i0,a)') solver, ': solout called ', ncalls, ' times'
    call check(err == 1, solver//' solout: status not 1')
    call check(ncalls >= 10, solver//' solout: too few calls')
    call check(numbers_consecutive, solver//' solout: step numbers not consecutive')
    call check(steps_join, solver//' solout: a step does not start where the previous one ended')
    call check(first_x == 0._R8, solver//' solout: first call not at the start time')
    call check(abs(last_x - TEND) <= 1.e-12_R8*TEND, solver//' solout: last call not at the end time')
    call check(all(last_y == y), solver//' solout: last solution differs from the one returned')
    call check(error_ratio(y, [cos(TEND), -sin(TEND)], [RT, RT], [AT, AT]) <= MAX_RATIO, &
               solver//' solout: error above the bound')
  end subroutine check_solout

  subroutine check_hmax(solver)
    character(len=*), intent(in) :: solver
    real(R8) :: y(2), t
    integer :: err
    call integrate_oscillator(solver, huge(1._R8), y, t, err, hmax_arg=HMAX)
    write(*,'(a10,a,es10.3,a,es10.3)') solver, ': longest step ', max_step, ' with hmax ', HMAX
    call check(err == 1, solver//' hmax: status not 1')
    call check(max_step <= HMAX*(1._R8 + 1.e-12_R8), solver//' hmax: a step longer than hmax')
    call check(ncalls - 1 >= nint(TEND/HMAX), solver//' hmax: fewer steps than TEND/hmax')
    call check(error_ratio(y, [cos(TEND), -sin(TEND)], [RT, RT], [AT, AT]) <= MAX_RATIO, &
               solver//' hmax: error above the bound')
  end subroutine check_hmax

  subroutine check_interrupt(solver)
    character(len=*), intent(in) :: solver
    real(R8) :: y(2), t
    integer :: err
    call integrate_oscillator(solver, STOP_TIME, y, t, err)
    write(*,'(a10,a,f8.4,a,i0)') solver, ': interrupted at t = ', t, ', status ', err
    call check(err == 2, solver//' interrupt: status not 2')
    call check(t >= STOP_TIME .and. t < TEND, solver//' interrupt: not stopped after the requested time')
    call check(last_x == t, solver//' interrupt: returned time differs from the last solout call')
    call check(error_ratio(y, [cos(t), -sin(t)], [RT, RT], [AT, AT]) <= MAX_RATIO, &
               solver//' interrupt: error above the bound')
  end subroutine check_interrupt

  !> 20 steps are too few for Robertson's problem up to 4e10 (the oscillator up to 1000 for H-dopri5)
  subroutine check_step_limit(solver)
    character(len=*), intent(in) :: solver
    real(R8) :: y(NROB), t, tout
    integer :: err
    if (solver == 'H-dopri5') then
      call setup_odesolver(N=2, solver=solver, RT=[RT, RT], AT=[AT, AT], IOPT=[20, 0, 0])
      y(1:2) = [1._R8, 0._R8]
      t = 0._R8
      tout = 1000._R8
      call run_odesolver(2, t, tout, y(1:2), oscillator, no_jacobian, 0, err)
    else
      call setup_odesolver(N=NROB, solver=solver, RT=spread(1.e-6_R8, 1, NROB), AT=[1.e-12_R8, 1.e-18_R8, 1.e-12_R8], &
                           IOPT=[20, 0, 0])
      y = [1._R8, 0._R8, 0._R8]
      t = 0._R8
      tout = ROB_TIMES(size(ROB_TIMES))
      call run_odesolver(NROB, t, tout, y, robertson, no_jacobian, 0, err)
    end if
    write(*,'(a10,a,i0)') solver, ': 20-step limit, status ', err
    call check(.not. solver_succeeded(solver, err), solver//' step limit: success reported with too few steps')
  end subroutine check_step_limit

end program options_test
