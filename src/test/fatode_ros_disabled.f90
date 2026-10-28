! In the builds that disable FATODE's Rosenbrock solvers (ifx with OpenMP, see CMakeLists.txt) setup_odesolver must
! refuse the solver named on the command line: it stops with an error and never reaches the last line. Built and run
! only in those builds, once for each of ros2, ros3, ros4, rodas3 and rodas4.

program fatode_ros_disabled_test
  use oslo
  implicit none
  character(len=10) :: solver

  call get_command_argument(1, solver)
  call setup_odesolver(N=3, solver=trim(solver))
  write(*,'(a)') 'setup_odesolver accepted '//trim(solver)

end program fatode_ros_disabled_test
