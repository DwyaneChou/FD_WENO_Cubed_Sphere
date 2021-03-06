module test_case_mod
  use parameters_mod
  use mesh_mod
  use ghost_mod
  use stat_mod
  use projection_mod
  use spatial_operators_mod
  use diag_mod
  implicit none
  
  contains
    
  subroutine initTestCase
    integer i,j,iPatch

    real, dimension(:,:,:), pointer :: phiG
    real, dimension(:,:,:), pointer :: u
    real, dimension(:,:,:), pointer :: v
    
    phiG => stat(0)%phiG
    u    => stat(0)%u   
    v    => stat(0)%v   
    
    if(case_num == 2)then
      print*,''
      print*,'test case 2 is selected'
      call case2(stat(0))
    elseif(case_num==5)then
      print*,''
      print*,'test case 5 is selected'
      call case5(stat(0))
    elseif(case_num == 6)then
      print*,''
      print*,'test case 6 is selected'
      call case6(stat(0))
    endif
    
    ! initialize phiG
    stat(0)%phiG = stat(0)%phi * mesh%sqrtG
    
    do iPatch = ifs, ife
      do j = jcs, jce
        do i = ics, ice
          call covProjSphere2Plane    (stat(0)%u (i,j,iPatch), stat(0)%v (i,j,iPatch), stat(0)%zonal_wind(i,j,iPatch), stat(0)%meridional_wind(i,j,iPatch), mesh%matrixIA(:,:,i,j,iPatch), mesh%matrixG(:,:,i,j,iPatch))
          call contravProjSphere2Plane(stat(0)%uC(i,j,iPatch), stat(0)%vC(i,j,iPatch), stat(0)%zonal_wind(i,j,iPatch), stat(0)%meridional_wind(i,j,iPatch), mesh%matrixIA(:,:,i,j,iPatch))
        enddo
      enddo
    enddo
      
    print*,''
    print*,'max/min value of phi : ',maxval(stat(0)%phi),minval(stat(0)%phi)
    print*,'max/min value of u   : ',maxval(stat(0)%u),minval(stat(0)%u)
    print*,'max/min value of v   : ',maxval(stat(0)%v),minval(stat(0)%v)
    print*,'max/min value of uc  : ',maxval(stat(0)%uc),minval(stat(0)%uc)
    print*,'max/min value of vc  : ',maxval(stat(0)%vc),minval(stat(0)%vc)
    print*,'max/min value of us  : ',maxval(stat(0)%zonal_wind),minval(stat(0)%zonal_wind)
    print*,'max/min value of vs  : ',maxval(stat(0)%meridional_wind),minval(stat(0)%meridional_wind)
  end subroutine initTestCase

  ! Global steady flow
  subroutine case2(stat)
    type(stat_field), intent(inout) :: stat
    real,dimension(ics:ice,jcs:jce,ifs:ife) :: u,v,phi                 ! working array
    
    real    :: u0
    real    :: gh0 = 29400.
    real    :: gh
    
    real,dimension(ics:ice,jcs:jce,ifs:ife) :: sinlat ! working array
    real,dimension(ics:ice,jcs:jce,ifs:ife) :: coslat ! working array
    
    integer :: i,j,iPatch
    
    sinlat = mesh%sinlat
    coslat = mesh%coslat
    
    u0 = 2. * pi * radius / (12. * 86400.)
    
    do iPatch = ifs, ife
      do j = jcs, jce
        do i = ics, ice
          phi(i,j,iPatch) = gh0 - (radius * Omega * u0 + u0**2 / 2.) * sinlat(i,j,iPatch)**2
          u  (i,j,iPatch) = u0 * coslat(i,j,iPatch)
        enddo
      enddo
    enddo
    
    v = 0.
    
    stat%zonal_wind      = u
    stat%meridional_wind = v
    stat%phi             = phi
    mesh%phis            = 0.
    
  end subroutine case2
  
  ! Isolated mountain
  subroutine case5(stat)
    type(stat_field), intent(inout) :: stat
    
    real,parameter :: hs0      = 2000.
    real,parameter :: u0       = 20.
    real,parameter :: alpha    = 0.
    real,parameter :: gh0      = 5960.*gravity
                                    
    real                                    :: rr
    real                                    :: labmda_c
    real                                    :: theta_c
    real,dimension(ics:ice,jcs:jce,ifs:ife) :: u,v,phi,ghs ! working array
    real,dimension(ics:ice,jcs:jce,ifs:ife) :: longitude  ! working array
    real,dimension(ics:ice,jcs:jce,ifs:ife) :: latitude   ! working array
    real,dimension(ics:ice,jcs:jce,ifs:ife) :: r
    
    integer :: i,j,iPatch ! working variable
    
    rr       = pi/9.
    labmda_c = 1.5*pi
    theta_c  = pi/6.
    
    longitude = mesh%lon
    latitude  = mesh%lat
    
    where(longitude<0) longitude = 2. * pi + longitude
    
    ghs = 0.
    
    do iPatch = ifs, ife
      do j = jcs, jce
          do i = ics, ice
              r  (i,j,iPatch) = sqrt(min(rr**2,(longitude(i,j,iPatch)-labmda_c)**2+(latitude(i,j,iPatch)-theta_c)**2))
              ghs(i,j,iPatch) = gravity*hs0*(1.-r(i,j,iPatch)/rr)
              u  (i,j,iPatch) = u0*(cos(latitude(i,j,iPatch))*cos(alpha)+cos(longitude(i,j,iPatch))*sin(latitude(i,j,iPatch))*sin(alpha))
              v  (i,j,iPatch) = -u0*sin(longitude(i,j,iPatch))*sin(alpha)
              phi(i,j,iPatch) = gh0 - (radius*Omega*u0 + u0**2/2.d0)*(-cos(longitude(i,j,iPatch))*cos(latitude(i,j,iPatch))*sin(alpha) + sin(latitude(i,j,iPatch))*cos(alpha))**2 - ghs(i,j,iPatch)
          enddo
      enddo
    enddo
    
    stat%zonal_wind      = u
    stat%meridional_wind = v
    stat%phi             = phi
    mesh%phis            = ghs
    
  end subroutine case5
  
  ! Rossby-Haurwitz wave with wavenumber 4
  subroutine case6(stat)
    type(stat_field), intent(inout) :: stat
    real,parameter              :: omg  = 7.848d-6         ! angular velocity of RH wave
    real,parameter              :: R    = 4              ! wave number of RH wave
    real,parameter              :: h0   = 8000.          ! wave number of RH wave
    
    real,dimension(ics:ice,jcs:jce,ifs:ife) :: u,v,phi                 ! working array
    real,dimension(ics:ice,jcs:jce,ifs:ife) :: u1,u2,u3                ! working array
    real,dimension(ics:ice,jcs:jce,ifs:ife) :: AA1,Ac,A21,A22,A23,Ah   ! working array
    real,dimension(ics:ice,jcs:jce,ifs:ife) :: Bc,BB1,BB2,Bh           ! working array
    real,dimension(ics:ice,jcs:jce,ifs:ife) :: CC,CC1,CC2,Ch           ! working array
    real,dimension(ics:ice,jcs:jce,ifs:ife) :: sinlat                  ! working array
    real,dimension(ics:ice,jcs:jce,ifs:ife) :: coslat                  ! working array
    real,dimension(ics:ice,jcs:jce,ifs:ife) :: longitude               ! working array
    
    integer                                 :: i,j,iPatch              ! working variable
    
    sinlat    = mesh%sinlat
    coslat    = mesh%coslat
    longitude = mesh%lon
    
    do iPatch = ifs, ife
      do j = jcs, jce
        do i = ics, ice
          u1(i,j,iPatch) = coslat(i,j,iPatch)
          u2(i,j,iPatch) = R*coslat(i,j,iPatch)**(R-1)*sinlat(i,j,iPatch)**2*cos(R*longitude(i,j,iPatch))
          u3(i,j,iPatch) = coslat(i,j,iPatch)**(R+1)*cos(R*longitude(i,j,iPatch))
          u (i,j,iPatch) = radius*omg*(u1(i,j,iPatch)+u2(i,j,iPatch)-u3(i,j,iPatch))
          
          v (i,j,iPatch) = -radius*omg*R*coslat(i,j,iPatch)**(R-1)*sinlat(i,j,iPatch)*sin(R*longitude(i,j,iPatch))
          
          AA1 (i,j,iPatch) = omg*0.5*(2.*Omega+omg)*coslat(i,j,iPatch)**2
          Ac  (i,j,iPatch) = 0.25*omg**2
          A21 (i,j,iPatch) = (R+1.)*coslat(i,j,iPatch)**(2.*R+2.)
          A22 (i,j,iPatch) = (2.*R**2-R-2.)*coslat(i,j,iPatch)**(2.*R)
          A23 (i,j,iPatch) = 2.*R**2*coslat(i,j,iPatch)**(2.*R-2)
          Ah  (i,j,iPatch) = AA1(i,j,iPatch)+Ac(i,j,iPatch)*(A21(i,j,iPatch)+A22(i,j,iPatch)-A23(i,j,iPatch))
          
          Bc  (i,j,iPatch) = 2.*(Omega+omg)*omg/((R+1)*(R+2))*coslat(i,j,iPatch)**R
          BB1 (i,j,iPatch) = R**2+2.*R+2.
          BB2 (i,j,iPatch) = (R+1.)**2.*coslat(i,j,iPatch)**2.
          Bh  (i,j,iPatch) = Bc(i,j,iPatch)*(BB1(i,j,iPatch)-BB2(i,j,iPatch))
          
          CC  (i,j,iPatch) = 0.25*omg**2*coslat(i,j,iPatch)**(2.*R)
          CC1 (i,j,iPatch) = (R+1.)*coslat(i,j,iPatch)**2;
          CC2 (i,j,iPatch) = R+2.
          Ch  (i,j,iPatch) = CC(i,j,iPatch)*(CC1(i,j,iPatch)-CC2(i,j,iPatch))
          
          phi (i,j,iPatch) = gravity*h0+radius**2*(Ah(i,j,iPatch) + Bh(i,j,iPatch)*cos(R*longitude(i,j,iPatch)) + Ch(i,j,iPatch)*cos(2.0*R*longitude(i,j,iPatch)))
        enddo
      enddo
    enddo
    
    stat%zonal_wind      = u
    stat%meridional_wind = v
    stat%phi             = phi
    mesh%phis            = 0.
  end subroutine case6
  
end module test_case_mod
    