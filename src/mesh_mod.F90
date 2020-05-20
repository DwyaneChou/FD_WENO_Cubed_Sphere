MODULE mesh_mod
  use constants_mod
  use parameters_mod
  use projection_mod
  
  implicit none
  
  ! coordinate
  type mesh_info
    real, dimension(:,:,:    ), allocatable :: xi       ! central angle on x direction for cells on each patch, unit: radian
    real, dimension(:,:,:    ), allocatable :: eta      ! central angle on y direction for cells on each patch, unit: radian
    real, dimension(:,:,:    ), allocatable :: x        ! x coordinate in cartesian
    real, dimension(:,:,:    ), allocatable :: y        ! y coordinate in cartesian
    real, dimension(:,:,:    ), allocatable :: z        ! z coordinate in cartesian
    real, dimension(:,:,:    ), allocatable :: lon      ! longitude on cells
    real, dimension(:,:,:    ), allocatable :: lat      ! latitude on cells
    real, dimension(:,:,:    ), allocatable :: sqrtG    ! jacobian of Transformation, sqrt(G)
    real, dimension(:,:,:,:,:), allocatable :: jab      ! jacobian matrix between cartesian and local panel
    real, dimension(:,:,:,:,:), allocatable :: matrixG  ! horizontal metric Tensor, which transform covariant vectors to contravariant vectors
    real, dimension(:,:,:,:,:), allocatable :: matrixIG ! horizontal metric Tensor, which transform contravariant vectors to covariant vectors
    real, dimension(:,:,:,:,:), allocatable :: matrixA  ! horizontal metric Tensor, which transform 
    real, dimension(:,:,:,:,:), allocatable :: matrixIA ! horizontal metric Tensor, which transform
    
    real, dimension(:,:,:    ), allocatable :: f       ! Coriolis parameter
    
    real, dimension(:,:,:    ), allocatable :: sinlon  ! sin(longitude)
    real, dimension(:,:,:    ), allocatable :: coslon  ! cos(longitude)
    real, dimension(:,:,:    ), allocatable :: sinlat  ! sin(latitude)
    real, dimension(:,:,:    ), allocatable :: coslat  ! cos(latitude)
    
    real, dimension(:,:,:    ), allocatable :: phis    ! surface geopotential height
    
    real, dimension(:,:,:    ), allocatable :: areaCell

    real, dimension(:,:,:    ), allocatable :: x_ext        ! Extended x coordinate in cartesian
    real, dimension(:,:,:    ), allocatable :: y_ext        ! Extended y coordinate in cartesian
    real, dimension(:,:,:    ), allocatable :: z_ext        ! Extended z coordinate in cartesian
    real, dimension(:,:,:    ), allocatable :: lon_ext      ! Extended longitude on cells
    real, dimension(:,:,:    ), allocatable :: lat_ext      ! Extended latitude on cells
    real, dimension(:,:,:    ), allocatable :: xi_ext       ! Extended xi on cells
    real, dimension(:,:,:    ), allocatable :: eta_ext      ! Extended eta on cells
    real, dimension(:,:,:    ), allocatable :: sqrtG_ext    ! Extended jacobian of Transformation, sqrt(G)
    real, dimension(:,:,:,:,:), allocatable :: jab_ext      ! Extended jacobian matrix between cartesian and local panel
    real, dimension(:,:,:,:,:), allocatable :: matrixG_ext  ! Extended horizontal metric Tensor, which transform covariant vectors to contravariant vectors
    real, dimension(:,:,:,:,:), allocatable :: matrixIG_ext ! Extended horizontal metric Tensor, which transform contravariant vectors to covariant vectors
    real, dimension(:,:,:,:,:), allocatable :: matrixA_ext  ! Extended horizontal metric Tensor, which transform 
    real, dimension(:,:,:,:,:), allocatable :: matrixIA_ext ! Extended horizontal metric Tensor, which transform
  end type mesh_info
  
  type(mesh_info), target :: mesh
  
  contains
  
  subroutine initMesh
    integer :: iPV, jPV, iCell, jCell, iPatch, iVar, iDOF, jDOF
    integer :: iPVs, iPVe, jPVs, jPVe
    
    ! Allocate arrays in structures
    allocate( mesh%x        (      ics:ice, jcs:jce, ifs:ife) )
    allocate( mesh%y        (      ics:ice, jcs:jce, ifs:ife) )
    allocate( mesh%z        (      ics:ice, jcs:jce, ifs:ife) )
    allocate( mesh%lon      (      ics:ice, jcs:jce, ifs:ife) )
    allocate( mesh%lat      (      ics:ice, jcs:jce, ifs:ife) )
    allocate( mesh%xi       (      ics:ice, jcs:jce, ifs:ife) )
    allocate( mesh%eta      (      ics:ice, jcs:jce, ifs:ife) )
    
    allocate( mesh%sqrtG    (      ics:ice, ics:ice, ifs:ife) )
    allocate( mesh%jab      (3, 2, ics:ice, ics:ice, ifs:ife) )
    allocate( mesh%matrixG  (2, 2, ics:ice, ics:ice, ifs:ife) )
    allocate( mesh%matrixIG (2, 2, ics:ice, ics:ice, ifs:ife) )
    allocate( mesh%matrixA  (2, 2, ics:ice, ics:ice, ifs:ife) )
    allocate( mesh%matrixIA (2, 2, ics:ice, ics:ice, ifs:ife) )
    
    allocate( mesh%f        (      ics:ice, ics:ice, ifs:ife) )
    
    allocate( mesh%sinlon   (      ics:ice, ics:ice, ifs:ife) )
    allocate( mesh%coslon   (      ics:ice, ics:ice, ifs:ife) )
    allocate( mesh%sinlat   (      ics:ice, ics:ice, ifs:ife) )
    allocate( mesh%coslat   (      ics:ice, ics:ice, ifs:ife) )
    
    allocate( mesh%phis     (      ics:ice, ics:ice, ifs:ife) )
    
    allocate( mesh%areaCell (      ics:ice, jcs:jce, ifs:ife) )
    
    allocate( mesh%x_ext        (      ids:ide, jds:jde, ifs:ife) )
    allocate( mesh%y_ext        (      ids:ide, jds:jde, ifs:ife) )
    allocate( mesh%z_ext        (      ids:ide, jds:jde, ifs:ife) )
    allocate( mesh%lon_ext      (      ids:ide, jds:jde, ifs:ife) )
    allocate( mesh%lat_ext      (      ids:ide, jds:jde, ifs:ife) )
    allocate( mesh%xi_ext       (      ids:ide, jds:jde, ifs:ife) )
    allocate( mesh%eta_ext      (      ids:ide, jds:jde, ifs:ife) )
    
    allocate( mesh%sqrtG_ext    (      ids:ide, ids:ide, ifs:ife) )
    allocate( mesh%jab_ext      (3, 2, ids:ide, ids:ide, ifs:ife) )
    allocate( mesh%matrixG_ext  (2, 2, ids:ide, ids:ide, ifs:ife) )
    allocate( mesh%matrixIG_ext (2, 2, ids:ide, ids:ide, ifs:ife) )
    allocate( mesh%matrixA_ext  (2, 2, ids:ide, ids:ide, ifs:ife) )
    allocate( mesh%matrixIA_ext (2, 2, ids:ide, ids:ide, ifs:ife) )
    
  end subroutine initMesh
END MODULE mesh_mod

