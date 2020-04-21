MODULE tend_mod
  use parameters_mod
  
  implicit none
  
  ! MCV basic definiton
  type tend_field
    real, dimension(:,:,:), allocatable :: phiG
    real, dimension(:,:,:), allocatable :: u   
    real, dimension(:,:,:), allocatable :: v   
  end type tend_field
  
  type(tend_field), dimension(:), target, allocatable :: tend ! allocated by n time points, which is used by temporal integration schemes
  
  contains
  
  subroutine initTend
    integer :: iT
    
    allocate( tend(-nIntegralSubSteps:1) )
    
    do iT = -nIntegralSubSteps, 1
      allocate(tend(iT)%phiG(ics:ice,jcs:jce,ifs:ife))
      allocate(tend(iT)%u   (ics:ice,jcs:jce,ifs:ife))
      allocate(tend(iT)%v   (ics:ice,jcs:jce,ifs:ife))
    enddo
    
  end subroutine initTend
  
END MODULE tend_mod

