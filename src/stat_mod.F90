MODULE stat_mod
  use parameters_mod
  
  implicit none
  
  ! MCV basic definiton
  type stat_field
    real, dimension(:,:,:), allocatable :: phiG  ! phi * sqrtG
    real, dimension(:,:,:), allocatable :: u     ! covariant u-wind
    real, dimension(:,:,:), allocatable :: v     ! covariant v-wind
    real, dimension(:,:,:), allocatable :: phi   ! geopotential height
    real, dimension(:,:,:), allocatable :: uC    ! contravariant u-wind
    real, dimension(:,:,:), allocatable :: vC    ! contravariant v-wind
    real, dimension(:,:,:), allocatable :: zonal_wind
    real, dimension(:,:,:), allocatable :: meridional_wind
    real, dimension(:,:,:), allocatable :: phit  ! phi + phis
  end type stat_field
  
  type(stat_field), dimension(:), target, allocatable :: stat  ! allocated by n time points, which is used by temporal integration schemes
  contains
  
  subroutine initStat
    integer :: iT
    
    allocate( stat(-nIntegralSubSteps:1) )
    
    do iT = -nIntegralSubSteps, 1
      allocate(stat(iT)%phiG           (ics:ice,jcs:jce,ifs:ife))
      allocate(stat(iT)%u              (ics:ice,jcs:jce,ifs:ife))
      allocate(stat(iT)%v              (ics:ice,jcs:jce,ifs:ife))
      allocate(stat(iT)%phi            (ics:ice,jcs:jce,ifs:ife))
      allocate(stat(iT)%uC             (ics:ice,jcs:jce,ifs:ife))
      allocate(stat(iT)%vC             (ics:ice,jcs:jce,ifs:ife))
      allocate(stat(iT)%zonal_wind     (ics:ice,jcs:jce,ifs:ife))
      allocate(stat(iT)%meridional_wind(ics:ice,jcs:jce,ifs:ife))
      allocate(stat(iT)%phit           (ics:ice,jcs:jce,ifs:ife))
      
      stat(iT)%phiG            = 0.
      stat(iT)%u               = 0.
      stat(iT)%v               = 0.
      stat(iT)%phi             = 0.
      stat(iT)%uC              = 0.
      stat(iT)%vC              = 0.
      stat(iT)%zonal_wind      = 0.
      stat(iT)%meridional_wind = 0.
      stat(iT)%phit            = 0.
    enddo
    
  end subroutine initStat
  
  subroutine copyStat(stat_out,stat_in)
    type(stat_field),intent(out) :: stat_out
    type(stat_field),intent(in ) :: stat_in
  
    stat_out%phiG            = stat_in%phiG           
    stat_out%u               = stat_in%u              
    stat_out%v               = stat_in%v              
    stat_out%phi             = stat_in%phi            
    stat_out%uC              = stat_in%uC     
    stat_out%vC              = stat_in%vC     
    stat_out%zonal_wind      = stat_in%zonal_wind     
    stat_out%meridional_wind = stat_in%meridional_wind
    stat_out%phit            = stat_in%phit
    
  end subroutine copyStat
  
END MODULE stat_mod

