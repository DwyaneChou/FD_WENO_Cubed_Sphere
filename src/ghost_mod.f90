MODULE ghost_mod
  use constants_mod
  use parameters_mod
  use projection_mod
  use mesh_mod
  use stat_mod
  
  implicit none
  
  contains
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
  
  SUBROUTINE CubedSphereFillJab(field)
  
    IMPLICIT NONE
  
    REAL, DIMENSION(3,2,ics:ice,jcs:jce,ifs:ife), INTENT(INOUT) :: field
    
    REAL, DIMENSION(3,2,its:ite,jts:jte,ifs:ife) :: field_inner
  
    ! Local variables
    INTEGER :: i
    
    !zarg = 0.0 !DBG
    field_inner = field(:,:,its:ite,jts:jte,ifs:ife)
  
    ! Equatorial panels
    !IF (np==1) THEN
       DO i=1,xhalo
          field(:,:,1-i    ,jts:jte,1) = field_inner(:,:,ite-i+1,jts:jte,4)  !exchange left
          field(:,:,ite+i  ,jts:jte,1) = field_inner(:,:,i      ,jts:jte,2)  !exchange right
          field(:,:,its:ite,1-i    ,1) = field_inner(:,:,its:ite,jte-i+1,6)  !exchange below
          field(:,:,its:ite,jte+i  ,1) = field_inner(:,:,its:ite,i      ,5)  !exchange over
       ENDDO
    !ELSE IF (np==2) THEN
       DO i=1,xhalo
          field(:,:,1-i    ,jts:jte,2) = field_inner(:,:,ite-i+1,jts:jte   ,1)  !exchange left
          field(:,:,ite+i  ,jts:jte,2) = field_inner(:,:,i      ,jts:jte   ,3)  !exchange right
          field(:,1,its:ite,1-i    ,2) =-field_inner(:,2,ite-i+1,jte:jts:-1,6)  !exchange below
          field(:,2,its:ite,1-i    ,2) = field_inner(:,1,ite-i+1,jte:jts:-1,6)  !exchange below
          field(:,1,its:ite,jte+i  ,2) = field_inner(:,2,ite-i+1,jts:jte   ,5)  !exchange over
          field(:,2,its:ite,jte+i  ,2) =-field_inner(:,1,ite-i+1,jts:jte   ,5)  !exchange over
       ENDDO
    !ELSE IF (np==3) THEN
       DO i=1,xhalo
          field(:,:,1-i    ,jts:jte,3) = field_inner(:,:,ite-i+1   ,jts:jte,2)  !exchange left
          field(:,:,ite+i  ,jts:jte,3) = field_inner(:,:,i         ,jts:jte,4)  !exchange right
          field(:,1,its:ite,1-i    ,3) =-field_inner(:,1,ite:its:-1,i      ,6)  !exchange below
          field(:,2,its:ite,1-i    ,3) =-field_inner(:,2,ite:its:-1,i      ,6)  !exchange below
          field(:,1,its:ite,jte+i  ,3) =-field_inner(:,1,ite:its:-1,jte-i+1,5)  !exchange over
          field(:,2,its:ite,jte+i  ,3) =-field_inner(:,2,ite:its:-1,jte-i+1,5)  !exchange over
       ENDDO
    !ELSE IF (np==4) THEN
       DO i=1,xhalo
          field(:,:,1-i    ,jts:jte,4) = field_inner(:,:,ite-i+1   ,jts:jte   ,3) !exchange left
          field(:,:,ite+i  ,jts:jte,4) = field_inner(:,:,i         ,jts:jte   ,1) !exchange right
          field(:,1,its:ite,1-i    ,4) = field_inner(:,2,i         ,jts:jte   ,6) !exchange below
          field(:,2,its:ite,1-i    ,4) =-field_inner(:,1,i         ,jts:jte   ,6) !exchange below
          field(:,1,its:ite,jte+i  ,4) =-field_inner(:,2,i         ,jte:jts:-1,5) !exchange over
          field(:,2,its:ite,jte+i  ,4) = field_inner(:,1,i         ,jte:jts:-1,5) !exchange over
       ENDDO
    ! Top panel
    !ELSE IF (np==5) THEN
       DO i=1,xhalo
          field(:,1,1-i    ,jts:jte,5) = field_inner(:,2,ite:its:-1,jte-i+1,4) !exchange left
          field(:,2,1-i    ,jts:jte,5) =-field_inner(:,1,ite:its:-1,jte-i+1,4) !exchange left
          field(:,1,ite+i  ,jts:jte,5) =-field_inner(:,2,its:ite   ,jte-i+1,2) !exchange right
          field(:,2,ite+i  ,jts:jte,5) = field_inner(:,1,its:ite   ,jte-i+1,2) !exchange right
          field(:,:,its:ite,1-i    ,5) = field_inner(:,:,its:ite   ,jte-i+1,1) !exchange below
          field(:,1,its:ite,jte+i  ,5) =-field_inner(:,1,ite:its:-1,jte-i+1,3) !exchange over
          field(:,2,its:ite,jte+i  ,5) =-field_inner(:,2,ite:its:-1,jte-i+1,3) !exchange over
       ENDDO
    ! Bottom panel
    !ELSE IF (np==6) THEN
       DO i=1,xhalo
          field(:,1,1-i    ,jts:jte,6) =-field_inner(:,2,its:ite   ,i,4) !exchange left
          field(:,2,1-i    ,jts:jte,6) = field_inner(:,1,its:ite   ,i,4) !exchange left
          field(:,1,ite+i  ,jts:jte,6) = field_inner(:,2,ite:its:-1,i,2) !exchange right
          field(:,2,ite+i  ,jts:jte,6) =-field_inner(:,1,ite:its:-1,i,2) !exchange right
          field(:,1,its:ite,1-i    ,6) =-field_inner(:,1,ite:its:-1,i,3) !exchange below
          field(:,2,its:ite,1-i    ,6) =-field_inner(:,2,ite:its:-1,i,3) !exchange below
          field(:,:,its:ite,jte+i  ,6) = field_inner(:,:,its:ite   ,i,1) !exchange over
       ENDDO
    !ELSE
    !   WRITE (*,*) 'Fatal error: In CubedSphereFillHalo'
    !   WRITE (*,*) 'Invalid panel id ', np
    !   STOP
    !ENDIF
  END SUBROUTINE CubedSphereFillJab
  
  ! Fill up halo with ghost points  
  subroutine fill_ghost(stat)
    type(stat_field), intent(inout) :: stat
    
    call CubedSphereFillHalo(stat%phiG           )
    call CubedSphereFillHalo(stat%zonal_wind     )
    call CubedSphereFillHalo(stat%meridional_wind)
    
  end subroutine fill_ghost
  
END MODULE ghost_mod