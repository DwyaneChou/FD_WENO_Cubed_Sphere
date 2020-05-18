MODULE ghost_mod
  use constants_mod
  use parameters_mod
  use projection_mod
  use mesh_mod
  use stat_mod
  
  implicit none
  
  ! ghost location
  type ghostLocation
    real   , dimension(:,:), allocatable :: X     ! X coordinate of ghost points
    real   , dimension(:,:), allocatable :: Y     ! Y coordinate of ghost points, not needed in interpolation
    real   , dimension(:,:), allocatable :: coef  ! Position in cell
    integer, dimension(:,:), allocatable :: iref  ! Index of reference point
  end type ghostLocation
  
  type(ghostLocation) :: ghost
  contains
  
  subroutine initGhost
    integer :: i,j
    integer :: iPatch ! computational patch
    integer :: gpatch ! ghost patch
    integer :: ihalo  ! halo index
    real    :: lambda, theta
    real    :: coef
    integer :: iref
    
    real    :: X_RAW(Nx+1)
    
    allocate(ghost%X   (its:ite,xhalo))
    allocate(ghost%Y   (jts:jte,xhalo))
    allocate(ghost%coef(its:ite,xhalo))
    allocate(ghost%iref(its:ite,xhalo))
    
    ! Calculate ghost point locations
    iPatch = 1
    gpatch = 5
    do j = jte+1, jce
      do i = its, ite
        ihalo = j - jte
        call pointProjPlane2Sphere(lambda          , theta           , mesh%x(i,j,iPatch), mesh%y(i,j,iPatch), iPatch)
        call pointProjSphere2Plane(ghost%X(i,ihalo), ghost%Y(i,ihalo), lambda            , theta             , gPatch)
      enddo
    enddo
    
    do i = 1, Nx
      X_RAW(i) = mesh%x(i,1,1)
    enddo
    
    ! Calculate reference points
    ! Now apply linear interpolation to obtain edge components
    DO j = 1, xhalo
      ! Reset the reference index
      iref = 1
      
      DO i = its, ite
        !
        ! Find reference points
        !
        DO WHILE (ghost%X(i,j) > X_RAW(iref))
          iref = iref + 1
        ENDDO
        
        IF ((ghost%X(i,j) >= X_RAW(iref-1)) .AND. (ghost%X(i,j) <= X_RAW(iref  )))THEN
            
          coef = ghost%X(i,j) - X_RAW(iref-1)
          
          ghost%coef(i,j) = coef
          ghost%iref(i,j) = iref - 1
        ELSE
          print*,'Ghost cell location error'
          stop
        ENDIF
      ENDDO
    ENDDO
    
  end subroutine initGhost
  
  !------------------------------------------------------------------------------
  ! SUBROUTINE CubedSphereFillHalo
  !
  ! Description:
  !   Recompute the cubed sphere data storage array, with the addition of a
  !   halo region around the specified panel.
  !------------------------------------------------------------------------------
  SUBROUTINE CubedSphereFillHalo(field)
  
    IMPLICIT NONE
  
    REAL, DIMENSION(ics:ice,jcs:jce,ifs:ife), INTENT(INOUT) :: field
    
    REAL, DIMENSION(its:ite,jts:jte,ifs:ife) :: field_inner
  
    ! Local variables
    INTEGER :: i
    
    !zarg = 0.0 !DBG
    field_inner = field(its:ite,jts:jte,ifs:ife)
  
    ! Equatorial panels
    !IF (np==1) THEN
       DO i=1,xhalo
          field(1-i    ,jts:jte,1) = field_inner(ite-i+1,jts:jte,4)  !exchange left
          field(ite+i  ,jts:jte,1) = field_inner(i      ,jts:jte,2)  !exchange right
          field(its:ite,1-i    ,1) = field_inner(its:ite,jte-i+1,6)  !exchange below
          field(its:ite,jte+i  ,1) = field_inner(its:ite,i      ,5)  !exchange over
       ENDDO
    !ELSE IF (np==2) THEN
       DO i=1,xhalo
          field(1-i    ,jts:jte,2) = field_inner(ite-i+1,jts:jte   ,1)  !exchange left
          field(ite+i  ,jts:jte,2) = field_inner(i      ,jts:jte   ,3)  !exchange right
          field(its:ite,1-i    ,2) = field_inner(ite-i+1,jte:jts:-1,6)  !exchange below
          field(its:ite,jte+i  ,2) = field_inner(ite-i+1,jts:jte   ,5)  !exchange over
       ENDDO
    !ELSE IF (np==3) THEN
       DO i=1,xhalo
          field(1-i    ,jts:jte,3) = field_inner(ite-i+1   ,jts:jte,2)  !exchange left
          field(ite+i  ,jts:jte,3) = field_inner(i         ,jts:jte,4)  !exchange right
          field(its:ite,1-i    ,3) = field_inner(ite:its:-1,i      ,6)  !exchange below
          field(its:ite,jte+i  ,3) = field_inner(ite:its:-1,jte-i+1,5)  !exchange over
       ENDDO
    !ELSE IF (np==4) THEN
       DO i=1,xhalo
          field(1-i    ,jts:jte,4) = field_inner(ite-i+1   ,jts:jte   ,3) !exchange left
          field(ite+i  ,jts:jte,4) = field_inner(i         ,jts:jte   ,1) !exchange right
          field(its:ite,1-i    ,4) = field_inner(i         ,jts:jte   ,6) !exchange below
          field(its:ite,jte+i  ,4) = field_inner(i         ,jte:jts:-1,5) !exchange over
       ENDDO
    ! Top panel
    !ELSE IF (np==5) THEN
       DO i=1,xhalo
          field(1-i    ,jts:jte,5) = field_inner(ite:its:-1,jte-i+1,4) !exchange left
          field(ite+i  ,jts:jte,5) = field_inner(its:ite   ,jte-i+1,2) !exchange right
          field(its:ite,1-i    ,5) = field_inner(its:ite   ,jte-i+1,1) !exchange below
          field(its:ite,jte+i  ,5) = field_inner(ite:its:-1,jte-i+1,3) !exchange over
       ENDDO
    ! Bottom panel
    !ELSE IF (np==6) THEN
       DO i=1,xhalo
          field(1-i    ,jts:jte,6) = field_inner(its:ite   ,i,4) !exchange left
          field(ite+i  ,jts:jte,6) = field_inner(ite:its:-1,i,2) !exchange right
          field(its:ite,1-i    ,6) = field_inner(ite:its:-1,i,3) !exchange below
          field(its:ite,jte+i  ,6) = field_inner(its:ite   ,i,1) !exchange over
       ENDDO
    !ELSE
    !   WRITE (*,*) 'Fatal error: In CubedSphereFillHalo'
    !   WRITE (*,*) 'Invalid panel id ', np
    !   STOP
    !ENDIF
  END SUBROUTINE CubedSphereFillHalo
  
  subroutine CubedSphereFillGhost(field)
    real, intent(inout) :: field(ics:ice,jcs:jce,ifs:ife)
    
    real ghost_target(ics:ice,jcs:jce,ifs:ife)
    
    integer i,j,iPatch
    
    ghost_target                          = FillValue
    ghost_target(its:ite,jts:jte,ifs:ife) = field(its:ite,jts:jte,ifs:ife)
    
    call CubedSphereFillHalo(field)
    
    do iPatch = ifs, ife
      do j = 1, xhalo
        ! Left boundary
        call linear_interp(ghost_target(its-j,jts:jte,iPatch),field(its-j,jts:jte,iPatch),ghost%iref(:,j),ghost%coef(:,j))
        ! Rigth boundary
        call linear_interp(ghost_target(ite+j,jts:jte,iPatch),field(ite+j,jts:jte,iPatch),ghost%iref(:,j),ghost%coef(:,j))
        ! Top boundary
        call linear_interp(ghost_target(its:ite,jte+j,iPatch),field(its:ite,jte+j,iPatch),ghost%iref(:,j),ghost%coef(:,j))
        ! Bottom boundary
        call linear_interp(ghost_target(its:ite,jts-j,iPatch),field(its:ite,jts-j,iPatch),ghost%iref(:,j),ghost%coef(:,j))
      enddo
    enddo
    
    field = ghost_target
  
  end subroutine CubedSphereFillGhost
  
  subroutine linear_interp(dest,src,iref,coef)
    real   , intent(out) :: dest(its:ite)
    real   , intent(in ) :: src (its:ite)
    integer, intent(in ) :: iref(its:ite)
    real   , intent(in ) :: coef(its:ite)
    
    integer i
    integer P1,P2,P3,P4,P5
    real    q1,q2,q3,q4,q5
    real    dh
    real    C
    
    dh = 4. * dx
    
    do i = its,ite
      if(iref(i)<=2)then
        P1 = 1
        P2 = 2
        P3 = 3
        P4 = 4
        P5 = 5
        C  = coef(i) + (iref(i)-1) * dx
      elseif(iref(i)>=ite-2)then
        P1 = ite - 4
        P2 = ite - 3
        P3 = ite - 2
        P4 = ite - 1
        P5 = ite
        C  = coef(i) + (4-ite+iref(i)) * dx
      else
        P1 = iref(i) - 2
        P2 = iref(i) - 1
        P3 = iref(i)
        P4 = iref(i) + 1
        P5 = iref(i) + 2
        C  = coef(i) + 2. * dx
      endif
      
      !P1 = iref(i)
      !P2 = iref(i)+1
      
      q1 = src(P1)
      q2 = src(P2)
      q3 = src(P3)
      q4 = src(P4)
      q5 = src(P5)
      
      dest(i) = ( 3. * dh**4 * q1 &
              +        dh**3 * (-25. * q1 + 48.  * q2 - 36.  * q3 + 16. * q4 -  3. * q5) * C     &
              +   2. * dh**2 * ( 35. * q1 - 104. * q2 + 114. * q3 - 56. * q4 + 11. * q5) * C**2  &
              -   16.* dh    * (  5. * q1 - 18.  * q2 + 24.  * q3 - 14. * q4 + 3.  * q5) * C**3  &
              +   32.*          (      q1 - 4.   * q2 + 6.   * q3 - 4.  * q4 +       q5) * C**4 )&
              /(3. * dh**4)
    
    enddo
  end subroutine linear_interp
  
  ! Fill up halo with ghost points  
  subroutine fill_ghost(stat)
    type(stat_field), intent(inout) :: stat
    
    call CubedSphereFillGhost(stat%phi            )
    call CubedSphereFillGhost(stat%zonal_wind     )
    call CubedSphereFillGhost(stat%meridional_wind)
    
  end subroutine fill_ghost
  
!!------------------------------------------------------------------------------
!! SUBROUTINE CubedSphereFillHalo
!!
!! Description:
!!   Recompute the cubed sphere data storage array, with the addition of a
!!   halo region around the specified panel.
!!
!! Parameters:
!!   parg - Current panel values
!!   zarg (OUT) - Calculated panel values with halo/ghost region
!!   np - Panel number
!!   ncube - Dimension of the cubed sphere (# of grid lines)
!!   nhalo - Number of halo/ghost elements around each panel
!!------------------------------------------------------------------------------
!  SUBROUTINE CubedSphereFillHalo(parg, zarg, np, ncube, nhalo)
!
!    IMPLICIT NONE
!
!    INTEGER , INTENT(IN) :: np, ncube,nhalo
!
!    REAL , DIMENSION(ncube-1, ncube-1, 6), INTENT(IN) :: parg
!
!    REAL ,                                            &
!         DIMENSION(1-nhalo:ncube+nhalo-1, 1-nhalo:ncube+nhalo-1, 6), &
!         INTENT(OUT) :: zarg
!
!
!    ! Local variables
!    INTEGER                :: jh
!    !zarg = 0.0 !DBG
!    zarg(1:ncube-1,1:ncube-1,np) = parg(1:ncube-1,1:ncube-1,np)
!
!    zarg(1-nhalo:0,1-nhalo:0,np) = 0.0
!    zarg(1-nhalo:0,ncube:ncube+nhalo-1,np) = 0.0
!    zarg(ncube:ncube+nhalo-1,1-nhalo:0,np) = 0.0
!    zarg(ncube:ncube+nhalo-1,ncube:ncube+nhalo-1,np) = 0.0
!
!    ! Equatorial panels
!    IF (np==1) THEN
!       DO jh=1,nhalo
!          zarg(1-jh      ,1:ncube-1 ,1) = parg(ncube-jh  ,1:ncube-1 ,4)  !exchange left
!          zarg(ncube+jh-1,1:ncube-1 ,1) = parg(jh        ,1:ncube-1 ,2)  !exchange right
!          zarg(1:ncube-1 ,1-jh      ,1) = parg(1:ncube-1 ,ncube-jh  ,6)  !exchange below
!          zarg(1:ncube-1 ,ncube+jh-1,1) = parg(1:ncube-1 ,jh        ,5)  !exchange over
!       ENDDO
!
!    ELSE IF (np==2) THEN
!       DO jh=1,nhalo
!          zarg(1-jh      ,1:ncube-1 ,2) = parg(ncube-jh,1:ncube-1   ,1)  !exchange left
!          zarg(ncube+jh-1,1:ncube-1 ,2) = parg(jh      ,1:ncube-1   ,3)  !exchange right
!          zarg(1:ncube-1 ,1-jh      ,2) = parg(ncube-jh,ncube-1:1:-1,6)  !exchange below
!          zarg(1:ncube-1 ,ncube+jh-1,2) = parg(ncube-jh,1:ncube-1   ,5)  !exchange over
!       ENDDO
!
!    ELSE IF (np==3) THEN
!       DO jh=1,nhalo
!          zarg(ncube+jh-1,1:ncube-1 ,3) = parg(jh          ,1:ncube-1,4)  !exchange right
!          zarg(1-jh      ,1:ncube-1 ,3) = parg(ncube-jh    ,1:ncube-1,2)  !exchange left
!          zarg(1:ncube-1 ,1-jh      ,3) = parg(ncube-1:1:-1,jh       ,6)  !exchange below
!          zarg(1:ncube-1 ,ncube+jh-1,3) = parg(ncube-1:1:-1,ncube-jh ,5)  !exchange over
!       ENDDO
!
!    ELSE IF (np==4) THEN
!       DO jh=1,nhalo
!          zarg(1-jh      ,1:ncube-1 ,4) = parg(ncube-jh,1:ncube-1   ,3) !exchange left
!          zarg(ncube+jh-1,1:ncube-1 ,4) = parg(jh      ,1:ncube-1   ,1) !exchange right
!          zarg(1:ncube-1 ,1-jh      ,4) = parg(jh      ,1:ncube-1   ,6) !exchange below
!          zarg(1:ncube-1 ,ncube+jh-1,4) = parg(jh      ,ncube-1:1:-1,5) !exchange over
!       ENDDO
!
!    ! Bottom panel
!    ELSE IF (np==5) THEN
!       DO jh=1,nhalo
!          zarg(1-jh      ,1:ncube-1 ,5) = parg(ncube-1:1:-1,ncube-jh,4) !exchange left
!          zarg(ncube+jh-1,1:ncube-1 ,5) = parg(1:ncube-1   ,ncube-jh,2) !exchange right
!          zarg(1:ncube-1 ,1-jh      ,5) = parg(1:ncube-1   ,ncube-jh,1) !exchange below
!          zarg(1:ncube-1 ,ncube+jh-1,5) = parg(ncube-1:1:-1,ncube-jh,3) !exchange over
!       ENDDO
!
!    ! Top panel
!    ELSE IF (np==6) THEN
!       DO jh=1,nhalo
!          zarg(1-jh      ,1:ncube-1 ,6) = parg(1:ncube-1   ,jh,4) !exchange left
!          zarg(ncube+jh-1,1:ncube-1 ,6) = parg(ncube-1:1:-1,jh,2) !exchange right
!          zarg(1:ncube-1 ,1-jh      ,6) = parg(ncube-1:1:-1,jh,3) !exchange below
!          zarg(1:ncube-1 ,ncube+jh-1,6) = parg(1:ncube-1   ,jh,1) !exchange over
!       ENDDO
!
!    ELSE
!       WRITE (*,*) 'Fatal error: In CubedSphereFillHalo'
!       WRITE (*,*) 'Invalid panel id ', np
!       STOP
!    ENDIF
!
!  END SUBROUTINE CubedSphereFillHalo

!------------------------------------------------------------------------------
! SUBROUTINE CubedSphereFillHalo_Linear
!
! Description:
!   Recompute the cubed sphere data storage array, with the addition of a
!   2-element halo region around the specified panel.  Use linear order
!   interpolation to translate between panels.
!
! Parameters:
!   parg - Current panel values
!   zarg (OUT) - Calculated panel values with halo/ghost region
!   np - Panel number
!   ncube - Dimension of the cubed sphere (# of grid lines)
!------------------------------------------------------------------------------
!  SUBROUTINE CubedSphereFillHalo_Linear(parg, zarg, np, ncube)
!
!!    USE CubedSphereTrans  ! Cubed sphere transforms
!
!    IMPLICIT NONE
!
!    INTEGER , PARAMETER :: nhalo = 2
!
!    INTEGER , INTENT(IN) :: np, ncube
!
!    REAL , DIMENSION(ncube-1, ncube-1, 6), INTENT(IN) :: parg
!
!    REAL ,                                            &
!         DIMENSION(1-nhalo:ncube+nhalo-1, 1-nhalo:ncube+nhalo-1, 6), &
!         INTENT(OUT) :: zarg
!
!
!    ! Local variables
!    INTEGER  :: ii, iref, jj, imin, imax
!    REAL     :: width, beta, a, newbeta
!
!    REAL    , DIMENSION(0:ncube, nhalo) :: prealpha
!    REAL    , DIMENSION(0:ncube, nhalo) :: newalpha
!
!    REAL , &
!         DIMENSION(1-nhalo:ncube+nhalo-1, 1-nhalo:ncube+nhalo-1, 6) :: yarg
!    
!    real, parameter :: one = 1., pih = 0.5 * pi, piq = 0.25 * pi
!
!    ! Use 0.0 order interpolation to begin
!    CALL CubedSphereFillHalo(parg, yarg, np, ncube, nhalo)
!
!    zarg(:,:,np) = yarg(:,:,np)
!
!    ! Calculate the overlapping alpha coordinates
!    width = pih / DBLE(ncube-1)
!
!    DO jj = 1, nhalo
!      DO ii = 0, ncube
!        prealpha(ii, jj) = width * (DBLE(ii-1) + 0.5) - piq
!        beta = - width * (DBLE(jj-1) + 0.5) - piq
!
!        CALL CubedSphereABPFromABP(prealpha(ii,jj), beta, 1, 5, &
!                                   newalpha(ii,jj), newbeta)
!      ENDDO
!    ENDDO
!
!    ! Now apply linear interpolation to obtain edge components
!    DO jj = 1, nhalo
!      ! Reset the reference index
!      iref = 2
!
!      ! Interpolation can be applied to more elements after first band
!      IF (jj == 1) THEN
!        imin = 1
!        imax = ncube-1
!      ELSE
!        imin = 0
!        imax = ncube
!      ENDIF
!
!      ! Apply linear interpolation
!      DO ii = imin, imax
!        DO WHILE ((iref .NE. ncube-1) .AND. &
!                  (newalpha(ii,jj) > prealpha(iref,jj)))
!          iref = iref + 1
!        ENDDO
!
!        IF ((newalpha(ii,jj) > prealpha(iref-1,jj)) .AND.    &
!            (newalpha(ii,jj) .LE. prealpha(iref  ,jj)))      &
!        THEN
!          a = (newalpha(ii,jj)   - prealpha(iref-1,jj)) / &
!              (prealpha(iref,jj) - prealpha(iref-1,jj))
!
!          IF ((a < 0.0) .OR. (a > one)) THEN
!            WRITE (*,*) 'FAIL in CubedSphereFillHalo_Linear'
!            WRITE (*,*) 'a out of bounds'
!            STOP
!          ENDIF
!
!          ! Bottom edge of panel
!          zarg(ii, 1-jj, np) =                   &
!            (one - a) * yarg(iref-1, 1-jj, np) + &
!                   a  * yarg(iref, 1-jj, np)
!
!          ! Left edge of panel
!          zarg(1-jj, ii, np) =                   &
!            (one - a) * yarg(1-jj, iref-1, np) + &
!                   a  * yarg(1-jj, iref, np)
!
!          ! Top edge of panel
!          zarg(ii, ncube+jj-1, np) =                   &
!            (one - a) * yarg(iref-1, ncube+jj-1, np) + &
!                   a  * yarg(iref, ncube+jj-1, np)
!
!          ! Right edge of panel
!          zarg(ncube+jj-1, ii, np) =                   &
!            (one - a) * yarg(ncube+jj-1, iref-1, np) + &
!                   a  * yarg(ncube+jj-1, iref, np)
!
!        ELSE
!          WRITE (*,*) 'FAIL in CubedSphereFillHalo_Linear'
!          WRITE (*,*) 'ii: ', ii, ' jj: ', jj
!          WRITE (*,*) 'newalpha: ', newalpha(ii,jj)
!          WRITE (*,*) 'prealpha: ', prealpha(iref-1,jj), '-', prealpha(iref,jj)
!          STOP
!        ENDIF
!      ENDDO
!    ENDDO
!
!    ! Fill in corner bits
!    zarg(0, 0, np) =                         &
!      0.25 * (zarg(1,0,np) + zarg(0,1,np) + &
!               zarg(-1,0,np) + zarg(0,-1,np))
!    zarg(0, ncube, np) =                                 &
!      0.25 * (zarg(0,ncube-1,np) + zarg(0,ncube+1,np) + &
!               zarg(-1,ncube,np)  + zarg(1,ncube,np))
!    zarg(ncube, 0, np) =                                 &
!      0.25 * (zarg(ncube-1,0,np) + zarg(ncube+1,0,np) + &
!               zarg(ncube,-1,np)  + zarg(ncube,1,np))
!    zarg(ncube, ncube, np) =                                     &
!      0.25 * (zarg(ncube-1,ncube,np) + zarg(ncube+1,ncube,np) + &
!               zarg(ncube,ncube-1,np) + zarg(ncube,ncube+1,np))
!
!  END SUBROUTINE CubedSphereFillHalo_Linear
!
!!--------------------------------------------------------------------------------
!! SUBROUTINE CubedSphereFillHalo_Linear_extended
!!
!! Same as CubedSphereFillHalo_Linear but it also fills the halo i<1 and i>ncube-1
!!
!!---------------------------------------------------------------------------------
!  SUBROUTINE CubedSphereFillHalo_Linear_extended(parg, parg_halo, np, ncube,nhalo)
!
!!    USE CubedSphereTrans  ! Cubed sphere transforms
!
!    IMPLICIT NONE
!
!    INTEGER , INTENT(IN) :: nhalo
!    INTEGER , INTENT(IN) :: np, ncube
!
!    REAL , DIMENSION(ncube-1, ncube-1, 6), INTENT(IN) :: parg
!
!    REAL ,                                            &
!         DIMENSION(1-nhalo:ncube+nhalo-1, 1-nhalo:ncube+nhalo-1), &
!         INTENT(OUT) :: parg_halo
!
!    REAL ,                                            &
!         DIMENSION(-nhalo:ncube+nhalo, -nhalo:ncube+nhalo):: zarg
!
!
!
!    ! Local variables
!    INTEGER  :: ii, iref, jj, imin, imax
!    REAL     :: width,beta, a, newbeta
!
!    REAL    , DIMENSION(-nhalo:ncube+nhalo, nhalo+1) :: prealpha !changed compared to non-extended
!    REAL    , DIMENSION(-nhalo:ncube+nhalo, nhalo+1) :: newalpha !changed compared to non-extended
!
!    REAL , &
!         DIMENSION(-nhalo:ncube+nhalo, -nhalo:ncube+nhalo, 6) :: yarg
!    
!    real, parameter :: one = 1., pih = 0.5 * pi, piq = 0.25 * pi
!
!    ! Use 0.0 order interpolation to begin
!    CALL CubedSphereFillHalo(parg, yarg, np, ncube, nhalo+1)
!
!    zarg(:,:) = yarg(:,:,np)
!
!    ! Calculate the overlapping alpha coordinates
!    width = pih / DBLE(ncube-1)
!
!    newalpha = -999999999999.0
!    prealpha = -999999999999.0
!
!    DO jj = 1, nhalo+1
!!      DO ii = 0, ncube
!      DO ii = 1-jj, ncube-1+jj !changed compared to non-extended
!        prealpha(ii, jj) = width * (DBLE(ii-1) + 0.5) - piq
!        beta = - width * (DBLE(jj-1) + 0.5) - piq
!
!        CALL CubedSphereABPFromABP(prealpha(ii,jj), beta, 1, 5, &
!                                   newalpha(ii,jj), newbeta)
!      ENDDO
!    ENDDO
!
!    ! Now apply linear interpolation to obtain edge components
!    DO jj = 1, nhalo+1
!      ! Reset the reference index
!      iref = 3-jj 
!
!      ! Interpolation can be applied to more elements after first band
!      !
!      imin = 2-jj       !changed compared to non-extended
!      imax = ncube-2+jj !changed compared to non-extended
!      !
!      ! Apply linear interpolation
!      !
!      DO ii = imin, imax
!        DO WHILE ((iref .NE. ncube-1) .AND. &
!                  (newalpha(ii,jj) > prealpha(iref,jj)))
!          iref = iref + 1
!        ENDDO
!
!        IF ((newalpha(ii,jj) > prealpha(iref-1,jj)) .AND.    &
!            (newalpha(ii,jj) .LE. prealpha(iref  ,jj)))      &
!        THEN
!          a = (newalpha(ii,jj)   - prealpha(iref-1,jj)) / &
!              (prealpha(iref,jj) - prealpha(iref-1,jj))
!
!          IF ((a < 0.0) .OR. (a > one)) THEN
!            WRITE (*,*) 'FAIL in CubedSphereFillHalo_Linear'
!            WRITE (*,*) 'a out of bounds'
!            STOP
!          ENDIF
!
!          ! Bottom edge of panel
!          zarg(ii, 1-jj) =                   &
!            (one - a) * yarg(iref-1, 1-jj, np) + &
!                   a  * yarg(iref, 1-jj, np)
!
!          ! Left edge of panel
!          zarg(1-jj, ii) =                   &
!            (one - a) * yarg(1-jj, iref-1, np) + &
!                   a  * yarg(1-jj, iref, np)
!
!          ! Top edge of panel
!          zarg(ii, ncube+jj-1) =                   &
!            (one - a) * yarg(iref-1, ncube+jj-1, np) + &
!                   a  * yarg(iref, ncube+jj-1, np)
!
!          ! Right edge of panel
!          zarg(ncube+jj-1, ii) =                   &
!            (one - a) * yarg(ncube+jj-1, iref-1, np) + &
!                   a  * yarg(ncube+jj-1, iref, np)
!
!        ELSE
!          WRITE (*,*) 'FAIL in CubedSphereFillHalo_Linear'
!          WRITE (*,*) 'ii: ', ii, ' jj: ', jj
!          WRITE (*,*) 'newalpha: ', newalpha(ii,jj)
!          WRITE (*,*) 'prealpha: ', prealpha(iref-1,jj), '-', prealpha(iref,jj)
!          STOP
!        ENDIF
!      ENDDO
!    ENDDO
!
!    ! Fill in corner bits
!    DO ii=0,nhalo-1
!      !
!      ! Diagonal lower left
!      !
!      zarg(0-ii, 0-ii) =                         &
!           0.25 * (zarg( 1-ii,0-ii) + zarg(0-ii,1-ii ) + &
!                   zarg(-1-ii,0-ii) + zarg(0-ii,-1-ii))
!      !
!      ! Diagonal upper left
!      !
!      zarg(0-ii, ncube+ii) =                                 &
!      0.25 * (zarg(0-ii,ncube-1+ii) + zarg(0-ii,ncube+1+ii) + &
!              zarg(-1-ii,ncube+ii ) + zarg(1-ii,ncube+ii  ))
!      !
!      ! Diagonal lower right
!      !
!      zarg(ncube+ii, 0-ii) =                                 &
!           0.25 * (zarg(ncube-1+ii, 0-ii) + zarg(ncube+1+ii,0-ii) + &
!                   zarg(ncube+ii  ,-1-ii) + zarg(ncube+ii  ,1-ii))
!      !
!      ! Diagonal upper right
!      !
!      zarg(ncube+ii, ncube+ii) =                                     &
!           0.25 * (zarg(ncube-1+ii,ncube+ii  ) + zarg(ncube+1+ii,ncube+ii  ) + &
!                   zarg(ncube+ii  ,ncube-1+ii) + zarg(ncube+ii  ,ncube+1+ii))
!    END DO
!
!    parg_halo = zarg(1-nhalo:ncube+nhalo-1, 1-nhalo:ncube+nhalo-1)
!  END SUBROUTINE CubedSphereFillHalo_Linear_extended
!
!
!!------------------------------------------------------------------------------
!! SUBROUTINE CubedSphereFillHalo_Cubic
!!
!! Description:
!!   Recompute the cubed sphere data storage array, with the addition of a
!!   2-element halo region around the specified panel.  Use higher order 
!!   interpolation to translate between panels.
!!
!! Parameters:
!!   parg - Current panel values
!!   zarg (OUT) - Calculated panel values with halo/ghost region
!!   np - Panel number
!!   ncube - Dimension of the cubed sphere (# of grid lines)
!!------------------------------------------------------------------------------
!  SUBROUTINE CubedSphereFillHalo_Cubic(parg, zarg, np, ncube)
!
!!    USE CubedSphereTrans  ! Cubed sphere transforms
!!    USE MathUtils         ! Has function for 1D cubic interpolation
!
!    IMPLICIT NONE
!
!    INTEGER , PARAMETER :: nhalo = 2
!
!    INTEGER , INTENT(IN) :: np, ncube
!
!    REAL , DIMENSION(ncube-1, ncube-1, 6), INTENT(IN) :: parg
!
!    REAL ,                                            &
!         DIMENSION(1-nhalo:ncube+nhalo-1, 1-nhalo:ncube+nhalo-1, 6), &
!         INTENT(OUT) :: zarg
!
!
!    ! Local variables
!    INTEGER  :: ii, iref, ibaseref, jj, imin, imax
!    REAL     :: width, beta, newbeta
!
!    REAL    , DIMENSION(0:ncube, nhalo) :: prealpha
!    REAL    , DIMENSION(0:ncube, nhalo) :: newalpha
!
!    REAL , &
!         DIMENSION(1-nhalo:ncube+nhalo-1, 1-nhalo:ncube+nhalo-1, 6) :: yarg
!    
!    real, parameter :: one = 1., pih = 0.5 * pi, piq = 0.25 * pi
!
!    ! Use 0.0 order interpolation to begin
!    CALL CubedSphereFillHalo(parg, yarg, np, ncube, nhalo)
!
!    zarg(:,:,np) = yarg(:,:,np)
!
!    ! Calculate the overlapping alpha coordinates
!    width = pih / DBLE(ncube-1)
!
!    DO jj = 1, nhalo
!      DO ii = 0, ncube
!        !
!        ! alpha,beta for the cell center (extending the panel)
!        !
!        prealpha(ii, jj) = width * (DBLE(ii-1) + 0.5) - piq
!        beta = - width * (DBLE(jj-1) + 0.5) - piq
!
!        CALL CubedSphereABPFromABP(prealpha(ii,jj), beta, 1, 5, &
!                                   newalpha(ii,jj), newbeta)
!      ENDDO
!    ENDDO
!
!    ! Now apply cubic interpolation to obtain edge components
!    DO jj = 1, nhalo
!      ! Reset the reference index, which gives the element in newalpha that
!      ! is closest to ii, looking towards larger values of alpha.
!      iref = 2 
!
!      ! Interpolation can be applied to more elements after first band
!!      IF (jj == 1) THEN
!!        imin = 1
!!        imax = ncube-1
!!      ELSE
!        imin = 0
!        imax = ncube
!!      ENDIF
!
!      ! Apply cubic interpolation
!      DO ii = imin, imax
!        DO WHILE ((iref .NE. ncube-1) .AND. &
!                  (newalpha(ii,jj) > prealpha(iref,jj)))
!          iref = iref + 1
!        ENDDO
!
!        ! Smallest index for cubic interpolation - apply special consideration
!        IF (iref == 2) THEN
!          ibaseref = iref-1
!
!        ! Largest index for cubic interpolation - apply special consideration
!        ELSEIF (iref == ncube-1) THEN
!          ibaseref = iref-3
!
!        ! Normal range
!        ELSE
!          ibaseref = iref-2
!        ENDIF
!
!        ! Bottom edge of panel
!        zarg(ii, 1-jj, np) =                                &
!          CUBIC_EQUISPACE_INTERP(                           &
!            width, newalpha(ii,jj) - prealpha(ibaseref,jj), &
!            yarg(ibaseref:ibaseref+3, 1-jj, np))
!
!        ! Left edge of panel
!        zarg(1-jj, ii, np) =                                &
!          CUBIC_EQUISPACE_INTERP(                           &
!            width, newalpha(ii,jj) - prealpha(ibaseref,jj), &
!            yarg(1-jj, ibaseref:ibaseref+3, np))
!
!        ! Top edge of panel
!        zarg(ii, ncube+jj-1, np) =                          &
!          CUBIC_EQUISPACE_INTERP(                           &
!            width, newalpha(ii,jj) - prealpha(ibaseref,jj), &
!            yarg(ibaseref:ibaseref+3, ncube+jj-1, np))
!
!        ! Right edge of panel
!        zarg(ncube+jj-1, ii, np) =                          &
!          CUBIC_EQUISPACE_INTERP(                           &
!            width, newalpha(ii,jj) - prealpha(ibaseref,jj), &
!            yarg(ncube+jj-1, ibaseref:ibaseref+3, np))
!
!      ENDDO
!    ENDDO
!
!    ! Fill in corner bits
!    zarg(0, 0, np) =                         &
!      0.25 * (zarg(1,0,np) + zarg(0,1,np) + &
!               zarg(-1,0,np) + zarg(0,-1,np))
!    zarg(0, ncube, np) =                                 &
!      0.25 * (zarg(0,ncube-1,np) + zarg(0,ncube+1,np) + &
!               zarg(-1,ncube,np)  + zarg(1,ncube,np))
!    zarg(ncube, 0, np) =                                 &
!      0.25 * (zarg(ncube-1,0,np) + zarg(ncube+1,0,np) + &
!               zarg(ncube,-1,np)  + zarg(ncube,1,np))
!    zarg(ncube, ncube, np) =                                     &
!      0.25 * (zarg(ncube-1,ncube,np) + zarg(ncube+1,ncube,np) + &
!               zarg(ncube,ncube-1,np) + zarg(ncube,ncube+1,np))
!
!  END SUBROUTINE CubedSphereFillHalo_Cubic
!
!!------------------------------------------------------------------------------
!! SUBROUTINE CubedSphereABPFromABP
!!
!! Description:
!!   Determine the (alpha,beta,idest) coordinate of a source point on
!!   panel isource.
!!
!! Parameters:
!!   alpha_in - Alpha coordinate in
!!   beta_in - Beta coordinate in
!!   isource - Source panel
!!   idest - Destination panel
!!   alpha_out (OUT) - Alpha coordinate out
!!   beta_out (OUT) - Beta coordiante out
!!------------------------------------------------------------------------------
!  SUBROUTINE CubedSphereABPFromABP(alpha_in,  beta_in, isource, idest, &
!                                   alpha_out, beta_out)
!
!    IMPLICIT NONE
!
!    REAL    , INTENT(IN)  :: alpha_in, beta_in
!    INTEGER , INTENT(IN)  :: isource, idest
!    REAL    , INTENT(OUT) :: alpha_out, beta_out
!
!    ! Local variables
!    REAL     :: a1, b1
!    REAL     :: xx, yy, zz
!    REAL     :: sx, sy, sz
!    
!    real, parameter :: one = 1., pih = 0.5 * pi, piq = 0.25 * pi
!
!    ! Convert to relative Cartesian coordinates
!    a1 = TAN(alpha_in)
!    b1 = TAN(beta_in)
!
!    sz = (one + a1 * a1 + b1 * b1)**(-0.5)
!    sx = sz * a1
!    sy = sz * b1
!
!    ! Convert to full Cartesian coordinates
!    IF (isource == 6) THEN
!      yy = sx; xx = -sy; zz = sz
!
!    ELSEIF (isource == 5) THEN
!      yy = sx; xx = sy; zz = -sz
!
!    ELSEIF (isource == 1) THEN
!      yy = sx; zz = sy; xx = sz
!
!    ELSEIF (isource == 3) THEN
!      yy = -sx; zz = sy; xx = -sz
!
!    ELSEIF (isource == 2) THEN
!      xx = -sx; zz = sy; yy = sz
!
!    ELSEIF (isource == 4) THEN
!      xx = sx; zz = sy; yy = -sz
!
!    ELSE
!      WRITE(*,*) 'Fatal Error: Source panel invalid in CubedSphereABPFromABP'
!      WRITE(*,*) 'panel = ', isource
!      STOP
!    ENDIF
!
!    ! Convert to relative Cartesian coordinates on destination panel
!    IF (idest == 6) THEN
!      sx = yy; sy = -xx; sz = zz
!
!    ELSEIF (idest == 5) THEN
!      sx = yy; sy = xx; sz = -zz
!
!    ELSEIF (idest == 1) THEN
!      sx = yy; sy = zz; sz = xx
!
!    ELSEIF (idest == 3) THEN
!      sx = -yy; sy = zz; sz = -xx
!
!    ELSEIF (idest == 2) THEN
!      sx = -xx; sy = zz; sz = yy
!
!    ELSEIF (idest == 4) THEN
!      sx = xx; sy = zz; sz = -yy
!
!    ELSE
!      WRITE(*,*) 'Fatal Error: Dest panel invalid in CubedSphereABPFromABP'
!      WRITE(*,*) 'panel = ', idest
!      STOP
!    ENDIF
!    IF (sz < 0) THEN
!      WRITE(*,*) 'Fatal Error: In CubedSphereABPFromABP'
!      WRITE(*,*) 'Invalid relative Z coordinate'
!      STOP
!    ENDIF
!
!    ! Use panel information to calculate (alpha, beta) coords
!    alpha_out = ATAN(sx / sz)
!    beta_out = ATAN(sy / sz)
!
!  END SUBROUTINE
!    
!!------------------------------------------------------------------------------
!! FUNCTION CUBIC_EQUISPACE_INTERP
!!
!! Description:
!!   Apply cubic interpolation on the specified array of values, where all
!!   points are equally spaced.
!!
!! Parameters:
!!   dx - Spacing of points
!!   x - X coordinate where interpolation is to be applied
!!   y - Array of 4 values = f(x + k * dx) where k = 0,1,2,3
!!------------------------------------------------------------------------------
!  FUNCTION CUBIC_EQUISPACE_INTERP(dx, x, y)
!    
!    IMPLICIT NONE
!    
!    REAL  :: CUBIC_EQUISPACE_INTERP
!    REAL  :: dx, x
!    REAL , DIMENSION(1:4) :: y
!    
!    CUBIC_EQUISPACE_INTERP =                                                   &
!         (-y(1) / (6.0 * dx**3)) * (x - dx) * (x - 2.0 * dx) * (x - 3.0 * dx) + &
!         ( y(2) / (2.0 * dx**3)) * (x) * (x - 2.0 * dx) * (x - 3.0 * dx) +      &
!         (-y(3) / (2.0 * dx**3)) * (x) * (x - dx) * (x - 3.0 * dx) +            &
!         ( y(4) / (6.0 * dx**3)) * (x) * (x - dx) * (x - 2.0 * dx)
!    
!  END FUNCTION CUBIC_EQUISPACE_INTERP
END MODULE ghost_mod