module parameters_mod
  use constants_mod
  implicit none
  
  integer, parameter :: nvar = 3
  
  ! Namelist parameters
  ! time_settings
  integer       :: run_days
  integer       :: run_hours
  integer       :: run_minutes
  integer       :: run_seconds
  real          :: dt               ! time step
  integer       :: history_interval ! output interval in seconds
  
  character*200 :: integral_scheme
  
  ! Case select
  integer :: case_num
  
  ! Domain
  character*200 :: mesh_file
  real          :: dx        !  grid-spacing in the x-direction
  real          :: dy        !  grid-spacing in the y-direction
  
  integer :: xhalo  !  halo number of x-diretion
  integer :: yhalo  !  halo number of y-diretion
  
  ! dynamic options
  character*200 :: reconstruct_scheme = 'WENO'
  
  ! Index parameter
  integer :: ids      ! The starting index in the x-direction (Physical domain)
  integer :: ide      ! The ending index in the x-direction  (Physical domain)
  integer :: jds      ! The starting index in the y-direction  (Physical domain)
  integer :: jde      ! The ending index in the y-direction  (Physical domain)
  
  integer :: ics      ! The starting index in the x-direction (Physical cell/element domain)
  integer :: ice      ! The ending index in the x-direction  (Physical cell/element domain)
  integer :: jcs      ! The starting index in the y-direction  (Physical cell/element domain)
  integer :: jce      ! The ending index in the y-direction  (Physical cell/element domain)
  
  integer :: its      ! The starting index in the x-direction of cell in physical domain
  integer :: ite      ! The ending index in the x-direction of cell in physical domain
  integer :: jts      ! The starting index in the y-direction of cell in physical domain
  integer :: jte      ! The ending index in the y-direction of cell in physical domain
  
  integer :: ips      ! The starting index in the x-direction (PV domain)
  integer :: ipe      ! The ending index in the x-direction  (PV domain)
  integer :: jps      ! The starting index in the y-direction  (PV domain)
  integer :: jpe      ! The ending index in the y-direction  (PV domain)
  
  integer :: ifs      ! The starting index of patch(face)
  integer :: ife      ! The ending index of patch(face)
                      
  integer :: Nx       ! Element numbers in the x-direction
  integer :: Ny       ! Element numbers in the y-direction
  
  integer :: Nx_halo  ! Element numbers in the x-direction with halo
  integer :: Ny_halo  ! Element numbers in the y-directionwith halo
  
  integer :: nPVx     ! Point-value numbers in the x-direction
  integer :: nPVy     ! Point-value numbers in the y-direction
  
  integer :: nPVx_halo! Point-value numbers in the x-direction
  integer :: nPVy_halo! Point-value numbers in the y-direction
  
  integer :: Nlambda  ! grid points in the lambda direction
  integer :: Ntheta   ! grid points in the theta direction
  
  integer :: nPVHalo
  
  integer :: nIntegralSubSteps ! number of integral substeps in temporal integration scheme
  integer :: nsteps            ! total integral steps
  
  integer, parameter :: Nf = 6           ! Number of cube faces
  
  real    :: x_min = -45.   !  start location of x-direction
  real    :: x_max =  45.   !  end location of x-direction
  real    :: y_min = -45.   !  start location of y-direction
  real    :: y_max =  45.   !  end location of y-direction
  
  ! Model run time control variables
  integer :: total_run_time   ! total run time for this model in seconds, this variable is determined by run_days, run_hours ...
  integer :: total_run_steps  ! total run steps for this model in seconds, this variable is determined by total_run_time and dt
  
  namelist /time_settings/ dt               ,&
                           run_days         ,&
                           run_hours        ,&
                           run_minutes      ,&
                           run_seconds      ,&
                           history_interval ,&
                           integral_scheme
  
  namelist /case_select/   case_num
  
  namelist /domain/        mesh_file   ,&
                           xhalo       ,&
                           yhalo
  namelist /dynamic_opt/   reconstruct_scheme
  
  contains
  
  subroutine readNamelist
    
    open(1, file = 'namelist.input',status='old')
    read(1, nml  = time_settings)
    read(1, nml  = case_select  )
    read(1, nml  = domain       )
    read(1, nml  = dynamic_opt  )
    close(1)
    
  end subroutine readNamelist
  
  subroutine initParameters
    
    ! Setting default values
    dt    = 300.
    dx    = 2.
    dy    = 2.
    xhalo = 1
    yhalo = 1
    
    run_days         = 1
    run_hours        = 0
    run_minutes      = 0
    run_seconds      = 0
    history_interval = 360
    integral_scheme  = 'RK4'
    
    ! Read namelist
    call readNamelist
    
    ! Calculate total run time in seconds
    total_run_time  = run_days * 86400 + run_hours * 3600 + run_minutes * 60 + run_seconds
    
    ! Calculate total run steps
    total_run_steps = ceiling(total_run_time/dt)
    
    ! Setting the number of substeps in temporal integration scheme
    if(trim(adjustl(integral_scheme)) == 'RK3_TVD')then
      nIntegralSubSteps = 3
    elseif(trim(adjustl(integral_scheme)) == 'RK4')then
      nIntegralSubSteps = 4
    else
      stop 'Unknown integral scheme, please select from RK3_TVD or RK4 ...'
    endif
    
    nsteps = total_run_time / dt
    
  end subroutine initParameters
  
end module parameters_mod
    