MODULE mesh_mod
  use constants_mod
  use parameters_mod
  use projection_mod
  
  implicit none
  
  ! coordinate
  type mesh_info
    real, dimension(:,:,:    ), allocatable :: x        ! central angle on x direction for cells on each patch, unit: radian
    real, dimension(:,:,:    ), allocatable :: y        ! central angle on y direction for cells on each patch, unit: radian
    real, dimension(:,:,:    ), allocatable :: lon      ! longitude on cells
    real, dimension(:,:,:    ), allocatable :: lat      ! latitude on cells
    real, dimension(:,:,:    ), allocatable :: sqrtG    ! jacobian of Transformation, sqrt(G)
    real, dimension(:,:,:,:,:), allocatable :: matrixG  ! horizontal metric Tensor, which transform covariant vectors to contravariant vectors
    real, dimension(:,:,:,:,:), allocatable :: matrixIG ! horizontal metric Tensor, which transform contravariant vectors to covariant vectors
    real, dimension(:,:,:,:,:), allocatable :: matrixA  ! horizontal metric Tensor, which transform 
    real, dimension(:,:,:,:,:), allocatable :: matrixIA ! horizontal metric Tensor, which transform
    
    real, dimension(:,:,:    ), allocatable :: f       ! Coriolis parameter
    
    real, dimension(:,:,:    ), allocatable :: sinlon  ! sin(longitude)
    real, dimension(:,:,:    ), allocatable :: coslon  ! cos(longitude)
    real, dimension(:,:,:    ), allocatable :: sinlat  ! sin(latitude)
    real, dimension(:,:,:    ), allocatable :: coslat  ! cos(latitude)
    
    real, dimension(:,:,:    ), allocatable :: sinx    ! trigonometric function
    real, dimension(:,:,:    ), allocatable :: cosx    ! trigonometric function
    real, dimension(:,:,:    ), allocatable :: tanx    ! trigonometric function
    real, dimension(:,:,:    ), allocatable :: cotx    ! trigonometric function
    real, dimension(:,:,:    ), allocatable :: secx    ! trigonometric function
    real, dimension(:,:,:    ), allocatable :: cscx    ! trigonometric function
    
    real, dimension(:,:,:    ), allocatable :: siny    ! trigonometric function
    real, dimension(:,:,:    ), allocatable :: cosy    ! trigonometric function
    real, dimension(:,:,:    ), allocatable :: tany    ! trigonometric function
    real, dimension(:,:,:    ), allocatable :: coty    ! trigonometric function
    real, dimension(:,:,:    ), allocatable :: secy    ! trigonometric function
    real, dimension(:,:,:    ), allocatable :: cscy    ! trigonometric function
    
    real, dimension(:,:,:    ), allocatable :: phis    ! surface geopotential height
    
    real, dimension(:,:,:    ), allocatable :: areaCell

    real, dimension(:,:,:    ), allocatable :: x_ext        ! Extended central angle on x direction for cells on each patch, unit: radian
    real, dimension(:,:,:    ), allocatable :: y_ext        ! Extended central angle on y direction for cells on each patch, unit: radian
    real, dimension(:,:,:    ), allocatable :: lon_ext      ! Extended longitude on cells
    real, dimension(:,:,:    ), allocatable :: lat_ext      ! Extended latitude on cells
    real, dimension(:,:,:    ), allocatable :: sqrtG_ext    ! Extended jacobian of Transformation, sqrt(G)
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
    
    real    :: areaCell_temp(Nx,Ny)
    
    ! Allocate arrays in structures
    allocate( mesh%x        (      ics:ice, jcs:jce, ifs:ife) )
    allocate( mesh%y        (      ics:ice, jcs:jce, ifs:ife) )
    allocate( mesh%lon      (      ics:ice, jcs:jce, ifs:ife) )
    allocate( mesh%lat      (      ics:ice, jcs:jce, ifs:ife) )
    
    allocate( mesh%sqrtG    (      ics:ice, ics:ice, ifs:ife) )
    allocate( mesh%matrixG  (2, 2, ics:ice, ics:ice, ifs:ife) )
    allocate( mesh%matrixIG (2, 2, ics:ice, ics:ice, ifs:ife) )
    allocate( mesh%matrixA  (2, 2, ics:ice, ics:ice, ifs:ife) )
    allocate( mesh%matrixIA (2, 2, ics:ice, ics:ice, ifs:ife) )
    
    allocate( mesh%f        (      ics:ice, ics:ice, ifs:ife) )
    
    allocate( mesh%sinlon   (      ics:ice, ics:ice, ifs:ife) )
    allocate( mesh%coslon   (      ics:ice, ics:ice, ifs:ife) )
    allocate( mesh%sinlat   (      ics:ice, ics:ice, ifs:ife) )
    allocate( mesh%coslat   (      ics:ice, ics:ice, ifs:ife) )
    
    allocate( mesh%sinx     (      ics:ice, ics:ice, ifs:ife) )
    allocate( mesh%cosx     (      ics:ice, ics:ice, ifs:ife) )
    allocate( mesh%tanx     (      ics:ice, ics:ice, ifs:ife) )
    allocate( mesh%cotx     (      ics:ice, ics:ice, ifs:ife) )
    allocate( mesh%secx     (      ics:ice, ics:ice, ifs:ife) )
    allocate( mesh%cscx     (      ics:ice, ics:ice, ifs:ife) )
    allocate( mesh%siny     (      ics:ice, ics:ice, ifs:ife) )
    allocate( mesh%cosy     (      ics:ice, ics:ice, ifs:ife) )
    allocate( mesh%tany     (      ics:ice, ics:ice, ifs:ife) )
    allocate( mesh%coty     (      ics:ice, ics:ice, ifs:ife) )
    allocate( mesh%secy     (      ics:ice, ics:ice, ifs:ife) )
    allocate( mesh%cscy     (      ics:ice, ics:ice, ifs:ife) )
    
    allocate( mesh%phis     (      ics:ice, ics:ice, ifs:ife) )
    
    allocate( mesh%areaCell (      ics:ice, jcs:jce, ifs:ife) )
    
    allocate( mesh%x_ext        (      ids:ide, jds:jde, ifs:ife) )
    allocate( mesh%y_ext        (      ids:ide, jds:jde, ifs:ife) )
    allocate( mesh%lon_ext      (      ids:ide, jds:jde, ifs:ife) )
    allocate( mesh%lat_ext      (      ids:ide, jds:jde, ifs:ife) )
    
    allocate( mesh%sqrtG_ext    (      ids:ide, ids:ide, ifs:ife) )
    allocate( mesh%matrixG_ext  (2, 2, ids:ide, ids:ide, ifs:ife) )
    allocate( mesh%matrixIG_ext (2, 2, ids:ide, ids:ide, ifs:ife) )
    allocate( mesh%matrixA_ext  (2, 2, ids:ide, ids:ide, ifs:ife) )
    allocate( mesh%matrixIA_ext (2, 2, ids:ide, ids:ide, ifs:ife) )
    
    ! Calculate mesh infomation on VIA
    do iPatch = ifs, ife
      do jCell = jcs, jce
        do iCell = ics, ice
          mesh%x(iCell, jCell, iPatch) = (iCell - 0.5) * dx + x_min
          mesh%y(iCell, jCell, iPatch) = (jCell - 0.5) * dy + y_min
          
          call pointProjPlane2Sphere(mesh%lon(iCell, jCell, iPatch), mesh%lat(iCell, jCell, iPatch), &
                                     mesh%x  (iCell, jCell, iPatch), mesh%y  (iCell, jCell, iPatch), iPatch)
          
          mesh%sinlon(iCell, jCell, iPatch) = sin(mesh%lon(iCell, jCell, iPatch))
          mesh%coslon(iCell, jCell, iPatch) = cos(mesh%lon(iCell, jCell, iPatch))
          
          mesh%sinlat(iCell, jCell, iPatch) = sin(mesh%lat(iCell, jCell, iPatch))
          mesh%coslat(iCell, jCell, iPatch) = cos(mesh%lat(iCell, jCell, iPatch))
          
          mesh%sinx(iCell, jCell, iPatch) = sin(mesh%x(iCell, jCell, iPatch))
          mesh%cosx(iCell, jCell, iPatch) = cos(mesh%x(iCell, jCell, iPatch))
          mesh%tanx(iCell, jCell, iPatch) = tan(mesh%x(iCell, jCell, iPatch))
          mesh%cotx(iCell, jCell, iPatch) = 1. / mesh%tanx(iCell, jCell, iPatch)
          mesh%secx(iCell, jCell, iPatch) = 1. / mesh%cosx(iCell, jCell, iPatch)
          mesh%cscx(iCell, jCell, iPatch) = 1. / mesh%sinx(iCell, jCell, iPatch)
          mesh%siny(iCell, jCell, iPatch) = sin(mesh%y(iCell, jCell, iPatch))
          mesh%cosy(iCell, jCell, iPatch) = cos(mesh%y(iCell, jCell, iPatch))
          mesh%tany(iCell, jCell, iPatch) = tan(mesh%y(iCell, jCell, iPatch))
          mesh%coty(iCell, jCell, iPatch) = 1. / mesh%tany(iCell, jCell, iPatch)
          mesh%secy(iCell, jCell, iPatch) = 1. / mesh%cosy(iCell, jCell, iPatch)
          mesh%cscy(iCell, jCell, iPatch) = 1. / mesh%siny(iCell, jCell, iPatch)
          
          call calc_matrixG (mesh%matrixG (:, :, iCell, jCell, iPatch), mesh%x  (iCell, jCell, iPatch), mesh%y  (iCell, jCell, iPatch))
          call calc_matrixIG(mesh%matrixIG(:, :, iCell, jCell, iPatch), mesh%x  (iCell, jCell, iPatch), mesh%y  (iCell, jCell, iPatch))
          call calc_matrixA (mesh%matrixA (:, :, iCell, jCell, iPatch), mesh%lon(iCell, jCell, iPatch), mesh%lat(iCell, jCell, iPatch), iPatch)
          call calc_matrixIA(mesh%matrixIA(:, :, iCell, jCell, iPatch), mesh%lon(iCell, jCell, iPatch), mesh%lat(iCell, jCell, iPatch), iPatch)
          call calc_Jacobian(mesh%sqrtG   (      iCell, jCell, iPatch), mesh%x  (iCell, jCell, iPatch), mesh%y  (iCell, jCell, iPatch))
          
          mesh%f(iCell, jCell, iPatch) = 2. * Omega * mesh%sinlat(iCell, jCell, iPatch)
        end do
      end do
    end do
    
    ! Calculate areaCell
    mesh%areaCell = 0.
    call EquiangularAllAreas(Nx, areaCell_temp)
    areaCell_temp = areaCell_temp * radius**2
    
    do iPatch = ifs, ife
      mesh%areaCell(1:Nx,1:Ny,iPatch) = areaCell_temp
    enddo
    
    ! Calculate extended mesh infomation on PV
    do iPatch = ifs, ife
      do jCell = jds, jde
        do iCell = ids, ide
          mesh%x_ext(iCell, jCell, iPatch) = 0.5 * (iCell - 1) * dx + x_min
          mesh%y_ext(iCell, jCell, iPatch) = 0.5 * (jCell - 1) * dy + y_min
          
          call pointProjPlane2Sphere(mesh%lon_ext(iCell, jCell, iPatch), mesh%lat_ext(iCell, jCell, iPatch), &
                                     mesh%x_ext  (iCell, jCell, iPatch), mesh%y_ext  (iCell, jCell, iPatch), iPatch)

          call calc_matrixG (mesh%matrixG_ext (:, :, iCell, jCell, iPatch), mesh%x_ext  (iCell, jCell, iPatch), mesh%y_ext  (iCell, jCell, iPatch))
          call calc_matrixIG(mesh%matrixIG_ext(:, :, iCell, jCell, iPatch), mesh%x_ext  (iCell, jCell, iPatch), mesh%y_ext  (iCell, jCell, iPatch))
          call calc_matrixA (mesh%matrixA_ext (:, :, iCell, jCell, iPatch), mesh%lon_ext(iCell, jCell, iPatch), mesh%lat_ext(iCell, jCell, iPatch), iPatch)
          call calc_matrixIA(mesh%matrixIA_ext(:, :, iCell, jCell, iPatch), mesh%lon_ext(iCell, jCell, iPatch), mesh%lat_ext(iCell, jCell, iPatch), iPatch)
          call calc_Jacobian(mesh%sqrtG_ext   (      iCell, jCell, iPatch), mesh%x_ext  (iCell, jCell, iPatch), mesh%y_ext  (iCell, jCell, iPatch))
        end do
      end do
    end do
    
    !do iCell = ids,ide
    !  print*,mesh%sqrtG_ext(iCell,ids,1)-mesh%sqrtG_ext(iCell,ide,6),mesh%sqrtG_ext(iCell,ide,1)
    !enddo
    
  end subroutine initMesh
  
  !------------------------------------------------------------------------------
  ! SUBROUTINE EquiangularAllAreas
  !
  ! Description:
  !   Compute the area of all cubed sphere grid cells, storing the results in
  !   a two dimensional array.
  !
  ! Parameters: 
  !   icube - Cell number (Nx or Ny) of the cubed sphere
  !   dA (OUT) - Output array containing the area of all cubed sphere grid cells
  !------------------------------------------------------------------------------
  SUBROUTINE EquiangularAllAreas(icube, dA)
    IMPLICIT NONE
    
    INTEGER,                         INTENT(IN)  :: icube
    REAL   , DIMENSION(icube,icube), INTENT(OUT) :: dA
    
    ! Local variables
    INTEGER                           :: k, k1, k2
    REAL                              :: a1, a2, a3, a4
    REAL , DIMENSION(icube+1,icube+1) :: ang
    REAL , DIMENSION(icube+1)         :: gp
    
    !#ifdef DBG 
    REAL    :: dbg !DBG
    !#endif
    
    ! Recall that we are using equi-angular spherical gridding
    !   Compute the angle between equiangular cubed sphere projection grid lines.
    DO k = 1, icube+1
      gp(k) = -0.25 * pi + (pi/DBLE(2*(icube))) * DBLE(k-1)
    ENDDO
    
    DO k2=1,icube+1
      DO k1=1,icube+1
        ang(k1,k2) = ACOS(-SIN(gp(k1)) * SIN(gp(k2)))
      ENDDO
    ENDDO
    
    DO k2=1,icube
      DO k1=1,icube
        a1 =      ang(k1  , k2  )
        a2 = pi - ang(k1+1, k2  )
        a3 = pi - ang(k1  , k2+1)
        a4 =      ang(k1+1, k2+1)      
        ! area = r*r*(-2*pi+sum(interior angles))
        DA(k1,k2) = -2.0*pi+a1+a2+a3+a4
      ENDDO
    ENDDO
    
    ! Only for debugging - test consistency
    dbg = 0.0
    DO k2=1,icube
      DO k1=1,icube
        dbg = dbg + DA(k1,k2)
      ENDDO
    ENDDO
    
    !print*,''
    !print*,'total area error     : ', dbg - 4. * pi / 6. !DBG
  END SUBROUTINE EquiangularAllAreas
END MODULE mesh_mod

