! OSLO - ODE SoLver toolbOx
!
! ODE solvers are set to work with the following specifications:
!   - full Jacobian, either evaluated by the solver (IJAC=0) or supplied by the
!     caller through the `jac` argument of run_odesolver (IJAC=1)
!   - the mass-matrix is the identity.
! All the parameters in the "setup" subroutine and within each solver procedure are set accordingly.
!
!***************************************************************************************************
!
! INTEGER INPUT
!
! Hairer solvers
!
! radau/radau5 -> 5
! 1) /=0 -> JACOBIAN MATRIX TO HESSENBERG FORM
! 2) max steps number
! 3) max number of newton iterations
! 4) choice of step initial value
! 8) gustafsson vs classic

! rodas -> 3
! 1) max steps number
! 2) choice of method
! 3) gustafsson vs classic

! sdirk4 -> 4
! 1) /=0 -> JACOBIAN MATRIX TO HESSENBERG FORM
! 2) max steps number
! 3) max number of newton iterations
! 4) coefficients

! FATODE solvers

! rk -> 6
! 3) method
! 4) max steps number
! 5) max number of newton iterations
! 6) choice of newton initial value
! 10) error estimation strategy
! 11) gustafsson vs classic

! ros -> 2
! 3) method
! 4) max steps number

! sdirk -> 4
! 3) method
! 4) max steps number
! 5) max number of newton iterations
! 6) choice of newton initial value

module interface_definitions
  implicit none

  interface
    subroutine solout_if(NR, XOLD, X, Y, N, IRTRN)
      use, intrinsic :: iso_fortran_env, only : I4 => int32, R8 => real64
      implicit none
      integer  :: NR, N
      real(R8) :: X, XOLD, Y(N)
      integer  :: IRTRN
    end subroutine solout_if
  end interface

  !> Analytical Jacobian of `fcn`. The signature is Hairer's: RADAU5, RODAS and
  !> SDIRK4 all declare JAC exactly this way, so one interface serves all three.
  !>   DFY(i,j) = d f(i) / d y(j), stored full (the wrappers set MLJAC = N).
  interface
    subroutine jac_if(N, X, Y, DFY, LDFY, RPAR, IPAR)
      use, intrinsic :: iso_fortran_env, only : I4 => int32, R8 => real64
      implicit none
      integer,  intent(in)  :: N, LDFY
      real(R8), intent(in)  :: X, Y(N)
      real(R8), intent(out) :: DFY(LDFY,N)
      real(R8), intent(in)  :: RPAR(*)
      integer,  intent(in)  :: IPAR(*)
    end subroutine jac_if
  end interface

end module interface_definitions

module oslo
  use, intrinsic :: iso_fortran_env, only : I4 => int32, R8 => real64
  implicit none
  private
  public :: setup_odesolver
  public :: run_odesolver
  public :: no_jacobian

  !> Concrete procedure pointing to one of the subroutine realizations
  procedure(solver_if), pointer :: run_odesolver

  !> `jac` is the analytical Jacobian of `fcn` and `IJAC` says whether to use it:
  !>   IJAC = 0 -> the solver builds the Jacobian by finite differences and never
  !>               touches `jac`; pass the `no_jacobian` stub exported above.
  !>   IJAC = 1 -> the solver calls `jac`.
  !> Both are mandatory and both are plain pass-through: no wrapper inspects them
  !> beyond handing them to the underlying integrator. Only the Hairer solvers
  !> (H-radau5, H-rodas, H-sdirk4) honour them -- see the individual wrappers.
  abstract interface
    subroutine solver_if(n, t1, t2, var, fcn, jac, IJAC, err, hmax, solout)
      use, intrinsic :: iso_fortran_env, only : I4 => int32, R8 => real64
      use interface_definitions, only : solout_if, jac_if
      implicit none
      integer, intent(in)            :: n
      real(R8), intent(inout)        :: t1, t2
      real(R8), intent(inout)        :: var(n)
      procedure(jac_if)              :: jac
      integer, intent(in)            :: IJAC
      integer, intent(out)           :: err
      real(R8), intent(in), optional :: hmax
      procedure(solout_if), optional :: solout
      external :: fcn
    end subroutine solver_if
  end interface

  !> Module data

  real(R8), allocatable :: RTOL(:)
  real(R8), allocatable :: ATOL(:)

  integer, parameter :: NIOPT=3
  integer, allocatable :: IWORK_global(:)
  integer :: ITOL
  integer :: LWORK, LIWORK
  integer :: LRCONT     ! see sdirk4.f
  integer :: NSMAX      ! see radau.f
  integer, parameter :: NNZERO=19 ! see FATODE

# if defined(SUNDIALS)
  !> Per-call context handed to CVODE as user_data. The RHS callback reads fcn
  !> and neq from here instead of through host association, so concurrent calls
  !> (one per OpenMP thread) never see each other's state.
  abstract interface
    subroutine cvode_fcn_if(n, t, y, f)
      import :: R8
      integer,  intent(in)  :: n
      real(R8), intent(in)  :: t
      real(R8), intent(in)  :: y(n)
      real(R8), intent(out) :: f(n)
    end subroutine cvode_fcn_if
  end interface
  type :: cvode_ctx_t
    procedure(cvode_fcn_if), pointer, nopass :: fcn => null()
    integer :: neq = 0
  end type cvode_ctx_t
# endif

contains

  subroutine setup_odesolver(N,solver,RT,AT,IOPT)
    implicit none
    integer, intent(in)            :: N
    character(len=*), intent(in)   :: solver
    real(R8), intent(in), optional :: RT(N), AT(N)
    integer, intent(in), optional  :: IOPT(NIOPT)
    integer :: Nsteps, Nnewt, icnt

    ! Tollerances are defined as arrays -> ITOL = 1 (Hairer)
    ITOL = 1

    Nsteps = 0; Nnewt = 0; icnt = 0
    if (present(IOPT)) then
      Nsteps = IOPT(1)
      Nnewt = IOPT(2)
      icnt = IOPT(3)
    endif

    if (allocated(IWORK_global)) deallocate(IWORK_global)

    if (allocated(RTOL)) deallocate(RTOL)
    if (allocated(ATOL)) deallocate(ATOL)
    nullify(run_odesolver)
    allocate(RTOL(1:N)); allocate(ATOL(1:N))
    ! Default values
    RTOL = 1.D-4; ATOL = 1.D-10
    if (present(RT)) RTOL = RT
    if (present(AT)) ATOL = AT

    select case(solver)

    case('H-dopri5')
      LWORK = 8*N + 5*N + 21
      LIWORK = N + 21
      allocate(IWORK_global(LIWORK))
      IWORK_global = 0
      IWORK_global(1) = Nsteps
      run_odesolver => wrap_dopri5

    case('H-radau5')
      LWORK=4*(n)*(n)+12*(n)+20
      LIWORK=3*n+20
      allocate(IWORK_global(LIWORK))
      IWORK_global = 0
      IWORK_global(2) = Nsteps
      IWORK_global(3) = Nnewt
      IWORK_global(8) = icnt
      run_odesolver => wrap_radau5

    case('radau2a','lobatto3c','gauss','radau1a')
      allocate(IWORK_global(NNZERO+1))
      IWORK_global = 0
      IWORK_global(4) = Nsteps
      IWORK_global(5) = Nnewt
      IWORK_global(11) = icnt
      IWORK_global(10) = 1
      if(solver=='radau2a')IWORK_global(3)=1
      if(solver=='lobatto3c')IWORK_global(3)=2
      if(solver=='gauss')IWORK_global(3)=3
      if(solver=='radau1a')IWORK_global(3)=4
      run_odesolver => wrap_rk_FATODE

    case('H-rodas')
      LWORK = 2*N*N+14*N+20
      LIWORK= N+20
      allocate(IWORK_global(LIWORK))
      IWORK_global = 0
      IWORK_global(1) = Nsteps
      IWORK_global(2) = Nnewt
      IWORK_global(3) = icnt
      run_odesolver => wrap_rodas

    case('ros2','ros3','ros4','rodas3','rodas4')
      allocate(IWORK_global(NNZERO+1))
      IWORK_global = 0
      IWORK_global(4) = Nsteps
      if(solver=='ros2')IWORK_global(3)=1
      if(solver=='ros3')IWORK_global(3)=2
      if(solver=='ros4')IWORK_global(3)=3
      if(solver=='rodas3')IWORK_global(3)=4
      if(solver=='rodas4')IWORK_global(3)=5
      run_odesolver => wrap_ros_FATODE

    case('H-sdirk4')
      LWORK=2*N*N+12*N+7
      LIWORK=2*N+4
      LRCONT=5*n+2
      allocate(IWORK_global(LIWORK))
      IWORK_global = 0
      IWORK_global(2) = Nsteps
      IWORK_global(3) = Nnewt
      run_odesolver => wrap_sdirk4

    case('sdirk2a','sdirk2b','sdirk3a','sdirk4a','sdirk4b')
      allocate(IWORK_global(NNZERO+1))
      IWORK_global = 0
      IWORK_global(4) = Nsteps
      IWORK_global(5) = Nnewt
      if(solver=='sdirk2a')IWORK_global(3)=1
      if(solver=='sdirk2b')IWORK_global(3)=2
      if(solver=='sdirk3a')IWORK_global(3)=3
      if(solver=='sdirk4a')IWORK_global(3)=4
      if(solver=='sdirk4b')IWORK_global(3)=5
      run_odesolver => wrap_sdirk_FATODE

#   if defined(INTEL)
    case('dodesol')
      LWORK = (7+2*n)*n
      LIWORK = 128
      allocate(IWORK_global(LIWORK))
      IWORK_global = 0
      run_odesolver => wrap_dodesol
#   endif

#   if defined(SUNDIALS)
    case('cvode')
      allocate(IWORK_global(1))
      IWORK_global = 10000
      run_odesolver => wrap_cvode
#   endif

    case default

      write(*,*) trim(solver), ' is not a valid integrator'
      stop

    end select

  end subroutine setup_odesolver


# if defined(INTEL)
!> `jac`/`IJAC` are inert here: using them means switching from dodesol_mk52lfn
!> to the analytic-Jacobian entry point dodesol_mk52lfa, a different call.
recursive subroutine wrap_dodesol(n,t1,t2,var,fcn,jac,IJAC,ierr,hmax,solout)
  use interface_definitions, only : solout_if, jac_if
  implicit none
  integer, intent(in)            :: n
  real(R8), intent(inout)        :: t1, t2
  real(R8), intent(inout)        :: var(n)
  procedure(jac_if)              :: jac
  integer, intent(in)            :: IJAC
  integer, intent(out)           :: ierr
  real(R8), intent(in), optional :: hmax
  external :: fcn
  procedure(solout_if), optional :: solout
  ! specific
  external :: dodesol_mk52lfn
  real(R8) :: h
  integer :: kd(n)
  real(R8) :: hm, ep, tr
  real(R8) :: dpar(LWORK)
  integer :: ipar(LIWORK)
  
  h = 1.d-15
  hm = 1.d-15
  ep = minval(RTOL)
  tr = minval(ATOL)
  ipar = IWORK_global

  call dodesol_mk52lfn(ipar,n,t1,t2,var,fcn,h,hm,ep,tr,dpar,kd,ierr)

end subroutine wrap_dodesol
# endif


recursive subroutine wrap_sdirk4(n,t1,t2,var,fcn,jac,IJAC,IDID,hmax,solout)
  use interface_definitions, only : solout_if, jac_if
  implicit none
  integer, intent(in)            :: n
  real(R8), intent(inout)        :: t1, t2
  real(R8), intent(inout)        :: var(n)
  procedure(jac_if)              :: jac
  integer, intent(in)            :: IJAC
  integer, intent(out)           :: IDID
  real(R8), intent(in), optional :: hmax
  procedure(solout_if), optional :: solout
  external :: fcn
  ! specific
  external :: SDIRK4
  integer :: MLJAC, MUJAC, MLMAS, MUMAS, IMAS
  real(R8) :: h
  real(R8) :: WORK(LWORK)
  integer :: IWORK(LIWORK)
  integer :: NN,NN2,NN3,NN4,IOUT
  real(R8) :: XOLD,HSOL,RCONT(LRCONT-2)
  integer :: NFCN,NJAC,NSTEP,NACCPT,NREJCT,NDEC,NSOL

  !> Seed the initial step from hmax when the caller supplies one. With h=0 SDIRK4
  !> picks its own 1e-6 default, which is far too large for stiff small-scale ODEs.
  if (present(hmax)) then; h = hmax * 1.0d-2
  else;                    h = 0.D0
  endif
  IMAS = 0
  MLJAC = n
  MLMAS = n
  MUJAC = 0
  MUMAS = 0
  ! IOUT=0: solout is never called
  ! IOUT=1: solout is available for output
  IOUT = 1

  IWORK = IWORK_global
  WORK = 0d0
  WORK(1)=1.1d-19
  if (present(hmax)) WORK(7) = hmax

  if (present(solout)) then
  call SDIRK4(n,fcn,t1,var,t2,H,                       &
                RTOL,ATOL,ITOL,                        &
                jac,IJAC,MLJAC,MUJAC,                  &
                dummy,IMAS,MLMAS,MUMAS,                &
                solout,IOUT,                           &
                WORK,LWORK,IWORK,LIWORK,LRCONT,IDID,   &
                NN,NN2,NN3,NN4,XOLD,HSOL,RCONT,        &
                NFCN,NJAC,NSTEP,NACCPT,NREJCT,NDEC,NSOL)
  else
  call SDIRK4(n,fcn,t1,var,t2,H,                       &
                RTOL,ATOL,ITOL,                        &
                jac,IJAC,MLJAC,MUJAC,                  &
                dummy,IMAS,MLMAS,MUMAS,                &
                dummy,IOUT,                           &
                WORK,LWORK,IWORK,LIWORK,LRCONT,IDID,   &
                NN,NN2,NN3,NN4,XOLD,HSOL,RCONT,        &
                NFCN,NJAC,NSTEP,NACCPT,NREJCT,NDEC,NSOL)
  endif

end subroutine wrap_sdirk4


recursive subroutine wrap_radau5(n,t1,t2,var,fcn,jac,IJAC,IDID,hmax,solout)
  use interface_definitions, only : solout_if, jac_if
  implicit none
  integer, intent(in)            :: n
  real(R8), intent(inout)        :: t1, t2
  real(R8), intent(inout)        :: var(n)
  procedure(jac_if)              :: jac
  integer, intent(in)            :: IJAC
  integer, intent(out)           :: IDID
  real(R8), intent(in), optional :: hmax
  procedure(solout_if), optional :: solout
  external :: fcn
  ! specific
  external :: RADAU5
  integer :: MLJAC, MUJAC, MLMAS, MUMAS, IMAS
  real(R8) :: h
  real(R8) :: WORK(LWORK)
  integer :: IWORK(LIWORK)
  real(R8) :: RPAR(1)
  integer :: IPAR(1)
  real(R8) :: RTOL_(n), ATOL_(n)
  integer :: NN,NN2,NN3,NN4,IOUT
  real(R8) :: XSOL,HSOL,C2M1,C1M1
  integer :: MLE,MUE,MBJAC,MBB,MDIAG,MDIFF,MBDIAG

  h = 0.D0
  IMAS = 0
  MLJAC = n
  MLMAS = n
  MUJAC = 0
  MUMAS = 0
  IOUT = 0

  IWORK = IWORK_global
  WORK = 0d0
  WORK(1)=1.1d-19
  if (present(hmax)) WORK(7) = hmax

  RPAR=0d0
  IPAR=0
  
  ! Reassignment is necessary since in radau5 tollerances are modified during the calculation
  RTOL_ = RTOL
  ATOL_ = ATOL

  call RADAU5(n,fcn,t1,var,t2,H,                        &
                RTOL_,ATOL_,ITOL,                       &
                jac,IJAC,MLJAC,MUJAC,                   &
                dummy,IMAS,MLMAS,MUMAS,                 &
                dummy,IOUT,                            &
                WORK,LWORK,IWORK,LIWORK,RPAR,IPAR,IDID, &
                NN,NN2,NN3,NN4,XSOL,HSOL,C2M1,C1M1,     &
                MLE,MUE,MBJAC,MBB,MDIAG,MDIFF,MBDIAG)

end subroutine wrap_radau5


recursive subroutine wrap_rodas(n,t1,t2,var,fcn,jac,IJAC,IDID,hmax,solout)
  use interface_definitions, only : solout_if, jac_if
  implicit none
  integer, intent(in)            :: n
  real(R8), intent(inout)        :: t1, t2
  real(R8), intent(inout)        :: var(n)
  procedure(jac_if)              :: jac
  integer, intent(in)            :: IJAC
  integer, intent(out)           :: IDID
  real(R8), intent(in), optional :: hmax
  procedure(solout_if), optional :: solout
  external :: fcn
  ! specific
  external :: RODAS
  integer :: MLJAC, MUJAC, MLMAS, MUMAS, IMAS
  real(R8) :: h
  real(R8) :: WORK(LWORK)
  integer :: IWORK(LIWORK)
  real(R8) :: RPAR(1)
  integer :: IPAR(1)
  integer :: IFCN, IDFX, IOUT

  h = 0.D0
  IFCN = 0
  IDFX = 0
  IMAS = 0
  MLJAC = n
  MLMAS = n
  MUJAC = 0
  MUMAS = 0
  IOUT = 0
  RPAR = 0d0
  IPAR = 0
  
  IWORK = IWORK_global
  WORK = 0d0
  WORK(1)=1.1d-19
  if (present(hmax)) WORK(2) = hmax

  call RODAS(n,fcn,IFCN,t1,var,t2,H,                   &
                RTOL,ATOL,ITOL,                        &
                jac,IJAC,MLJAC,MUJAC,dummy,IDFX,       &
                dummy,IMAS,MLMAS,MUMAS,                &
                dummy,IOUT,                           &
                WORK,LWORK,IWORK,LIWORK,RPAR,IPAR,IDID)

end subroutine wrap_rodas


!> `jac`/`IJAC` are inert here: DOPRI5 is explicit and never forms a Jacobian.
recursive subroutine wrap_dopri5(n,t1,t2,var,fcn,jac,IJAC,IDID,hmax,solout)
  use interface_definitions, only : solout_if, jac_if
  implicit none
  integer, intent(in)            :: n
  real(R8), intent(inout)        :: t1, t2
  real(R8), intent(inout)        :: var(n)
  procedure(jac_if)              :: jac
  integer, intent(in)            :: IJAC
  integer, intent(out)           :: IDID
  real(R8), intent(in), optional :: hmax
  procedure(solout_if), optional :: solout
  external :: fcn
  external :: DOPRI5
  real(R8) :: WORK(LWORK)
  integer :: IWORK(LIWORK), IOUT
  real(R8) :: RTOL_(n), ATOL_(n)

  IWORK = IWORK_global
  WORK = 0d0
  if (present(hmax)) WORK(6) = hmax
  RTOL_ = RTOL; ATOL_ = ATOL
  IOUT = 1

  if (present(solout)) then
    call DOPRI5(n,fcn,t1,var,t2,RTOL_,ATOL_,ITOL, &
                solout,IOUT, &
                WORK,LWORK,IWORK,LIWORK,IDID)
  else
    call DOPRI5(n,fcn,t1,var,t2,RTOL_,ATOL_,ITOL, &
                dummy,IOUT, &
                WORK,LWORK,IWORK,LIWORK,IDID)
  end if

end subroutine wrap_dopri5


!> `jac`/`IJAC` are inert here: FATODE declares its Jacobian callback as
!> JAC(T,Y,Jac0), which is not Hairer's signature, so it cannot be forwarded.
recursive subroutine wrap_sdirk_FATODE(n,t1,t2,var,fcn,jac,IJAC,err,hmax,solout)
  use SDIRK_f90_Integrator
  use interface_definitions, only : solout_if, jac_if
  implicit none
  integer, intent(in)            :: n
  real(R8), intent(inout)        :: t1, t2
  real(R8), intent(inout)        :: var(n)
  procedure(jac_if)              :: jac
  integer, intent(in)            :: IJAC
  integer, intent(out)           :: err
  real(R8), intent(in), optional :: hmax
  procedure(solout_if), optional :: solout
  external :: fcn
  ! specific
  real(R8) :: RCNTRL(NNZERO+1), RSTATUS(NNZERO+1)
  integer  :: ICNTRL(NNZERO+1), ISTATUS(NNZERO+1)
   
  ICNTRL  = IWORK_global
  RCNTRL  = 0.0

  if (present(hmax)) RCNTRL(2) = hmax
  
  call SDIRK(N,NNZERO,T1,T2,VAR,RTOL,ATOL,fcn,dummy,  &
             RCNTRL,ICNTRL,RSTATUS,ISTATUS,err)

end subroutine wrap_sdirk_FATODE


!> `jac`/`IJAC` are inert here -- see wrap_sdirk_FATODE.
recursive subroutine wrap_ros_FATODE(n,t1,t2,var,fcn,jac,IJAC,err,hmax,solout)
  use ROS_f90_Integrator, only: Rosenbrock
  use interface_definitions, only : solout_if, jac_if
  implicit none
  integer, intent(in)            :: n
  real(R8), intent(inout)        :: t1, t2
  real(R8), intent(inout)        :: var(n)
  procedure(jac_if)              :: jac
  integer, intent(in)            :: IJAC
  integer, intent(out)           :: err
  real(R8), intent(in), optional :: hmax
  procedure(solout_if), optional :: solout
  external :: fcn
  ! specific
  real(R8) :: RCNTRL(NNZERO+1), RSTATUS(NNZERO+1)
  integer  :: ICNTRL(NNZERO+1), ISTATUS(NNZERO+1)
   
  ICNTRL  = IWORK_global
  RCNTRL  = 0.0
  RSTATUS = 0.0
  ISTATUS = 0

  if (present(hmax)) RCNTRL(2) = hmax

  call Rosenbrock(N,NNZERO,VAR,T1,T2,   &
          ATOL,RTOL, fcn, dummy,        &
          RCNTRL,ICNTRL,RSTATUS,ISTATUS,err)

end subroutine wrap_ros_FATODE


!> `jac`/`IJAC` are inert here -- see wrap_sdirk_FATODE.
recursive subroutine wrap_rk_FATODE(n,t1,t2,var,fcn,jac,IJAC,err,hmax,solout)
  use RK_f90_Integrator
  use interface_definitions, only : solout_if, jac_if
  implicit none
  integer, intent(in)            :: n
  real(R8), intent(inout)        :: t1, t2
  real(R8), intent(inout)        :: var(n)
  procedure(jac_if)              :: jac
  integer, intent(in)            :: IJAC
  integer, intent(out)           :: err
  real(R8), intent(in), optional :: hmax
  procedure(solout_if), optional :: solout
  external :: fcn
  ! specific
  real(R8) :: RCNTRL(NNZERO+1), RSTATUS(NNZERO+1)
  integer  :: ICNTRL(NNZERO+1), ISTATUS(NNZERO+1)
   
  ICNTRL  = IWORK_global
  RCNTRL  = 0.0
  if (present(hmax)) RCNTRL(2) = hmax

  call RungeKutta(  N, NNZERO, T1, T2, VAR, RTOL, ATOL, fcn, dummy, &
                    RCNTRL,ICNTRL,RSTATUS,ISTATUS,err )

end subroutine wrap_rk_FATODE

# if defined(SUNDIALS)
!> `jac`/`IJAC` are inert here: CVODES takes its Jacobian through FCVodeSetJacFn
!> as a C-interoperable callback, not through the argument list.
recursive subroutine wrap_cvode(n,t1,t2,var,fcn,jac,IJAC,err,hmax,solout)
  use, intrinsic :: iso_c_binding
  use fsundials_core_mod
  use fcvodes_mod                   ! Fortran interface to CVODES
  use fnvector_serial_mod           ! Fortran interface to serial N_Vector
  use fsunmatrix_dense_mod          ! Fortran interface to dense SUNMatrix
  use fsunlinsol_dense_mod          ! Fortran interface to dense SUNLinearSolver
  use fsunnonlinsol_newton_mod
  use interface_definitions, only : solout_if, jac_if
  implicit none
  integer, intent(in)            :: n
  real(R8), intent(inout)        :: t1, t2
  real(R8), intent(inout)        :: var(n)
  procedure(jac_if)              :: jac
  integer, intent(in)            :: IJAC
  integer, intent(out)           :: err
  real(R8), intent(in), optional :: hmax
  procedure(solout_if), optional :: solout
  external :: fcn
  ! specific
  type(N_Vector), pointer :: sunvec_y      ! sundials solution vector
  type(N_Vector), pointer :: sunvec_f      ! sundials solution vector
  type(N_Vector), pointer :: sunvec_av     ! sundials tolerance vector
  type(SUNMatrix), pointer :: sunmat_A      ! sundials matrix
  type(SUNLinearSolver), pointer :: sunlinsol_LS  ! sundials linear solver
  type(SUNNonLinearSolver), pointer :: sunnonlin_NLS ! sundials nonlinear solver
  type(c_ptr)                       :: cvode_mem     ! CVode memory
  type(c_ptr)                       :: sunctx        ! SUNDIALS simulation context
  type(cvode_ctx_t), target         :: ctx           ! RHS context, see cvode_rhs
  real(R8)                          :: fval(n), tret(1)
  integer(c_int)                    :: retval
  integer(c_int64_t)                :: neq

  neq = int(n, c_int64_t)

  retval = FSUNContext_Create(SUN_COMM_NULL, sunctx)

  sunvec_y => FN_VMake_Serial(neq, var, sunctx)
  sunvec_f => FN_VMake_Serial(neq, fval, sunctx)
  sunvec_av => FN_VMake_Serial(neq, ATOL, sunctx)

  ! Call FCVodeCreate and FCVodeInit to create and initialize CVode memory
  ctx%fcn => fcn
  ctx%neq = n

  cvode_mem = FCVodeCreate(CV_BDF, sunctx)
  retval = FCVodeInit(cvode_mem, c_funloc(cvode_rhs), t1, sunvec_y)
  retval = FCVodeSetUserData(cvode_mem, c_loc(ctx))

  ! Set tolerances and maxstep
  retval = FCVodeSVtolerances(cvode_mem, RTOL(1), sunvec_av)
  retval = FCVodeSetMaxNumSteps(cvode_mem, int(IWORK_global(1), c_long))

  ! Create dense SUNMatrix for use in linear solves
  sunmat_A => FSUNDenseMatrix(neq, neq, sunctx)
  ! Create dense SUNLinearSolver object
  sunlinsol_LS => FSUNLinSol_Dense(sunvec_y, sunmat_A, sunctx)
  ! Attach the matrix and linear solver
  retval = FCVodeSetLinearSolver(cvode_mem, sunlinsol_LS, sunmat_A)
  ! Use Newton non-linear solver
  sunnonlin_NLS => FSUNNonlinSol_Newton(sunvec_y, sunctx)
  retval = FCVodeSetNonlinearSolver(cvode_mem, sunnonlin_NLS)
  ! Assign Jacobian
  !retval = FCVodeSetJacFn(cvode_mem, c_funloc(jac))

  if (present(hmax)) &
    retval =  FCVodeSetMaxStep(cvode_mem, hmax)

  ! A STOP here would kill the whole run from inside the caller's parallel
  ! region: hand the CVODE status back instead (negative = failure).
  retval = FCVode(cvode_mem, t2, sunvec_y, tret(1), CV_NORMAL)
  err = int(retval)

  !call PrintFinalStats(cvode_mem)

  ! free memory (this runs once per call, so anything left here leaks per cell)
  call FCVodeFree(cvode_mem)
  retval = FSUNNonlinSolFree(sunnonlin_NLS)
  retval = FSUNLinSolFree(sunlinsol_LS)
  call FSUNMatDestroy(sunmat_A)
  call FN_VDestroy(sunvec_y)
  call FN_VDestroy(sunvec_f)
  call FN_VDestroy(sunvec_av)
  retval = FSUNContext_Free(sunctx)


contains


  ! ----------------------------------------------------------------
  ! PrintFinalStats
  !
  ! Print CVode statstics to standard out
  ! ----------------------------------------------------------------
  subroutine PrintFinalStats(cvode_mem)

    !======= Inclusions ===========
    use iso_c_binding
    use fcvodes_mod

    !======= Declarations =========
    implicit none

    type(c_ptr), intent(in) :: cvode_mem ! solver memory structure

    integer(c_int)  :: retval          ! error flag

    integer(c_long) :: nsteps(1)     ! num steps
    integer(c_long) :: nfe(1)        ! num function evals
    integer(c_long) :: netfails(1)   ! num error test fails
    integer(c_long) :: nniters(1)    ! nonlinear solver iterations
    integer(c_long) :: nncfails(1)   ! nonlinear solver fails
    integer(c_long) :: njacevals(1)  ! number of Jacobian evaluations
    integer(c_long) :: nluevals(1)   ! number of LU evals
    integer(c_long) :: ngevals(1)    ! number of root evals

    !======= Internals ============

    retval = FCVodeGetNumSteps(cvode_mem, nsteps)
    if (retval /= 0) then
      print *, 'Error in FCVodeGetNumSteps, retval = ', retval, '; halting'
      stop 1
    end if

    retval = FCVodeGetNumRhsEvals(cvode_mem, nfe)
    if (retval /= 0) then
      print *, 'Error in FCVodeGetNumRhsEvals, retval = ', retval, '; halting'
      stop 1
    end if

    retval = FCVodeGetNumLinSolvSetups(cvode_mem, nluevals)
    if (retval /= 0) then
      print *, 'Error in FCVodeGetNumLinSolvSetups, retval = ', retval, '; halting'
      stop 1
    end if

    retval = FCVodeGetNumErrTestFails(cvode_mem, netfails)
    if (retval /= 0) then
      print *, 'Error in FCVodeGetNumErrTestFails, retval = ', retval, '; halting'
      stop 1
    end if

    retval = FCVodeGetNumNonlinSolvIters(cvode_mem, nniters)
    if (retval /= 0) then
      print *, 'Error in FCVodeGetNumNonlinSolvIters, retval = ', retval, '; halting'
      stop 1
    end if

    retval = FCVodeGetNumNonlinSolvConvFails(cvode_mem, nncfails)
    if (retval /= 0) then
      print *, 'Error in FCVodeGetNumNonlinSolvConvFails, retval = ', retval, '; halting'
      stop 1
    end if

    retval = FCVodeGetNumJacEvals(cvode_mem, njacevals)
    if (retval /= 0) then
      print *, 'Error in FCVodeGetNumJacEvals, retval = ', retval, '; halting'
      stop 1
    end if

    retval = FCVodeGetNumGEvals(cvode_mem, ngevals)
    if (retval /= 0) then
      print *, 'Error in FCVodeGetNumGEvals, retval = ', retval, '; halting'
      stop 1
    end if

    print *, ' '
    print *, ' General Solver Stats:'
    print '(4x,A,i9)', 'Total internal steps taken    =', nsteps
    print '(4x,A,i9)', 'Total rhs function calls      =', nfe
    print '(4x,A,i9)', 'Total Jacobian function calls =', njacevals
    print '(4x,A,i9)', 'Total LU function calls       =', nluevals
    print '(4x,A,i9)', 'Num error test failures       =', netfails
    print '(4x,A,i9)', 'Num nonlinear solver iters    =', nniters
    print '(4x,A,i9)', 'Num nonlinear solver fails    =', nncfails
    print *, ' '

    return

  end subroutine PrintFinalStats

end subroutine wrap_cvode


!> CVODE RHS callback. It is a module procedure (not internal to wrap_cvode) and
!> takes everything it needs from user_data, so it is reentrant: a C function
!> pointer to an internal procedure carries its host link through a trampoline,
!> which breaks when several threads integrate at once.
recursive integer(c_int) function cvode_rhs(t, sunvec_y, sunvec_f, user_data) result(ierr) bind(C)
  use, intrinsic :: iso_c_binding
  use fsundials_core_mod
  use fnvector_serial_mod
  implicit none
  real(c_double), value :: t         ! current time
  type(N_Vector)        :: sunvec_y  ! solution N_Vector
  type(N_Vector)        :: sunvec_f  ! function N_Vector
  type(c_ptr), value    :: user_data ! cvode_ctx_t of the calling wrap_cvode
  type(cvode_ctx_t), pointer :: ctx
  real(c_double), pointer :: yval(:), fval(:)
  real(R8) :: tloc

  call c_f_pointer(user_data, ctx)
  yval => FN_VGetArrayPointer(sunvec_y)
  fval => FN_VGetArrayPointer(sunvec_f)

  tloc = t
  call ctx%fcn(ctx%neq, tloc, yval, fval)
  ! fcn flags an out-of-range state by returning -1 everywhere: recoverable error
  ierr = 0
  if (all(fval == -1.0d0)) ierr = 1

end function cvode_rhs
# endif


subroutine dummy()
endsubroutine


!> Stub to pass as `jac` when IJAC = 0, so that argument can stay mandatory.
!> The integrators never call it in that case; zeroing DFY keeps it harmless
!> if some future solver does.
subroutine no_jacobian(N, X, Y, DFY, LDFY, RPAR, IPAR)
  implicit none
  integer,  intent(in)  :: N, LDFY
  real(R8), intent(in)  :: X, Y(N)
  real(R8), intent(out) :: DFY(LDFY,N)
  real(R8), intent(in)  :: RPAR(*)
  integer,  intent(in)  :: IPAR(*)
  DFY = 0.0_R8
end subroutine no_jacobian


end module oslo
