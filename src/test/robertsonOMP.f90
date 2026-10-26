! Thread safety of every solver: NCOPIES integrations run concurrently from an OpenMP parallel region, as a solver of
! many cells calls OSLO. The copies end at different times (the cases), so that threads work on different states at
! the same moment: shared work arrays or COMMON blocks would mix them, which identical copies could hide. Every copy
! must return a successful status and the same solution, bit for bit, as a serial run of its case.
! The stiff solvers integrate Robertson's problem (and the Hairer ones once more with the analytical Jacobian), H-dopri5
! the oscillator. The test fails if it runs on one thread only.

program robertsonOMP_test
  use, intrinsic :: iso_c_binding, only : c_loc, c_intptr_t
  use omp_lib
  use oslo
  use oslo_test_support
  implicit none
  integer, parameter :: NCOPIES = 240
  real(R8), parameter :: RT(NROB) = 1.e-6_R8
  real(R8), parameter :: AT(NROB) = [1.e-12_R8, 1.e-18_R8, 1.e-12_R8]
  real(R8), parameter :: OSC_TIMES(6) = [1._R8, 2._R8, 3._R8, 5._R8, 8._R8, 13._R8]
  character(len=10), allocatable :: solvers(:)
  integer :: i

  write(*,'(a,i0,a)') 'Running on ', omp_get_max_threads(), ' threads'
  call check(omp_get_max_threads() > 1, 'more than one OpenMP thread (set OMP_NUM_THREADS)')

  solvers = stiff_solvers()
  do i = 1, size(solvers)
    call run_concurrent(trim(solvers(i)), 0)
  end do
  solvers = hairer_jacobian_solvers()
  do i = 1, size(solvers)
    call run_concurrent(trim(solvers(i)), 1)
  end do
  call run_concurrent('H-dopri5', 0)

  call test_summary()

contains

  !> Integrates case c from 0 to its end time: Robertson's problem up to ROB_TIMES(c), or the oscillator up to
  !> OSC_TIMES(c) for H-dopri5
  subroutine integrate_case(solver, ijac, c, y, err)
    character(len=*), intent(in) :: solver
    integer, intent(in) :: ijac, c
    real(R8), intent(out) :: y(NROB)
    integer, intent(out) :: err
    real(R8) :: t, tout
    t = 0._R8
    if (solver == 'H-dopri5') then
      y = [1._R8, 0._R8, 0._R8]
      tout = OSC_TIMES(c)
      call run_odesolver(2, t, tout, y(1:2), oscillator, no_jacobian, 0, err)
    else
      y = [1._R8, 0._R8, 0._R8]
      tout = ROB_TIMES(c)
      if (ijac == 1) then
        call run_odesolver(NROB, t, tout, y, robertson, robertson_jac, 1, err)
      else
        call run_odesolver(NROB, t, tout, y, robertson, no_jacobian, 0, err)
      end if
    end if
  end subroutine integrate_case

  subroutine run_concurrent(solver, ijac)
    character(len=*), intent(in) :: solver
    integer, intent(in) :: ijac
    integer, parameter :: NCASES = size(ROB_TIMES)
    real(R8) :: yserial(NROB,NCASES)
    real(R8), target :: y(NROB,NCOPIES)
    integer :: err(NCOPIES), thread(NCOPIES), err_serial, c, i, nbad_status, nbad_value
    character(len=64) :: label

    write(label,'(a,a,i0,a)') solver, ' (ijac=', ijac, ')'
    if (solver == 'H-dopri5') then
      call setup_odesolver(N=2, solver=solver, RT=RT(1:2), AT=[AT(1), AT(1)])
    else
      call setup_odesolver(N=NROB, solver=solver, RT=RT, AT=AT)
    end if

    do c = 1, NCASES
      call integrate_case(solver, ijac, c, yserial(:,c), err_serial)
      call check(solver_succeeded(solver, err_serial), trim(label)//': serial run unsuccessful')
    end do

    !$omp parallel do schedule(dynamic,1) default(shared) private(i)
    do i = 1, NCOPIES
      call integrate_case(solver, ijac, case_of(i), y(:,i), err(i))
      thread(i) = omp_get_thread_num()
    end do
    !$omp end parallel do

    nbad_status = 0
    nbad_value = 0
    do i = 1, NCOPIES
      if (.not. solver_succeeded(solver, err(i))) nbad_status = nbad_status + 1
      if (any(y(:,i) /= yserial(:,case_of(i)))) nbad_value = nbad_value + 1
    end do
    write(*,'(a24,a,i0,a,i0,a,i0,a)') trim(label), ': ', NCOPIES, ' copies, ', nbad_status, ' unsuccessful, ', &
                                       nbad_value, ' differing from the serial run'
    call check(nbad_status == 0, trim(label)//': concurrent runs unsuccessful')
    call check(nbad_value == 0, trim(label)//': concurrent runs differ from the serial run')
    if (nbad_value > 0) call explain_differences(solver, ijac, y, yserial, thread)
  end subroutine run_concurrent

  !> What the differing copies have in common: their case, their thread, the alignment of their solution vector.
  !> The copies are also run again without threads, in an array laid out as y: if these differ as well, the results
  !> depend on where the solution vector lies in memory, not on the threads
  subroutine explain_differences(solver, ijac, y, yserial, thread)
    character(len=*), intent(in) :: solver
    integer, intent(in) :: ijac
    real(R8), intent(in), target :: y(:,:)
    real(R8), intent(in) :: yserial(:,:)
    integer, intent(in) :: thread(:)
    real(R8), target :: z(NROB,NCOPIES)
    integer :: by_case(size(yserial,2)), by_thread(0:size(thread)-1), all_by_thread(0:size(thread)-1), by_align(0:7)
    integer :: i, c, nbad_layout, align, err
    real(R8) :: worst, serial_error, parallel_error
    by_case = 0; by_thread = 0; all_by_thread = 0; by_align = 0
    worst = 0._R8
    serial_error = 0._R8
    parallel_error = 0._R8
    do i = 1, NCOPIES
      all_by_thread(thread(i)) = all_by_thread(thread(i)) + 1
    end do
    do i = 1, NCOPIES
      if (all(y(:,i) == yserial(:,case_of(i)))) cycle
      align = int(mod(transfer(c_loc(y(1,i)), 0_c_intptr_t), 64_c_intptr_t)) / 8
      by_case(case_of(i)) = by_case(case_of(i)) + 1
      by_thread(thread(i)) = by_thread(thread(i)) + 1
      by_align(align) = by_align(align) + 1
      worst = max(worst, error_ratio(y(:,i), yserial(:,case_of(i)), RT, AT))
      if (solver /= 'H-dopri5') parallel_error = max(parallel_error, error_ratio(y(:,i), ROB_REF(:,case_of(i)), RT, AT))
    end do
    if (solver /= 'H-dopri5') then
      do c = 1, size(yserial,2)
        serial_error = max(serial_error, error_ratio(yserial(:,c), ROB_REF(:,c), RT, AT))
      end do
    end if
    nbad_layout = 0
    do i = 1, NCOPIES
      call integrate_case(solver, ijac, case_of(i), z(:,i), err)
      if (any(z(:,i) /= yserial(:,case_of(i)))) nbad_layout = nbad_layout + 1
    end do
    write(*,'(4x,a,es10.3)') 'largest difference, in units of the tolerance: ', worst
    write(*,'(4x,a,*(1x,i0))') 'differing copies by case:              ', by_case
    write(*,'(4x,a,*(1x,i0))') 'differing copies by thread:            ', by_thread(0:maxval(thread))
    write(*,'(4x,a,*(1x,i0))') 'all copies by thread:                  ', all_by_thread(0:maxval(thread))
    if (solver /= 'H-dopri5') then
      write(*,'(4x,a,es10.3)') 'error of the serial runs against the reference, in tolerances:      ', serial_error
      write(*,'(4x,a,es10.3)') 'error of the differing copies against the reference, in tolerances: ', parallel_error
    end if
    write(*,'(4x,a,*(1x,i0))') 'differing copies by address mod 64 / 8:', by_align
    write(*,'(4x,a,i0,a,i0)')  'copies run again without threads, differing: ', nbad_layout, ' of ', NCOPIES
  end subroutine explain_differences

  integer function case_of(i)
    integer, intent(in) :: i
    case_of = mod(i - 1, size(ROB_TIMES)) + 1
  end function case_of

end program robertsonOMP_test
