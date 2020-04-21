MODULE spatial_operators_mod
  use constants_mod
  use mesh_mod
  use parameters_mod
  use stat_mod
  use tend_mod
  use projection_mod
  
  implicit none
  
    contains
    
    subroutine spatial_operator(stat,tend)
      type(stat_field), target, intent(in   ) :: stat
      type(tend_field), target, intent(inout) :: tend
      
      real, dimension(ics:ice,jcs:jce,ifs:ife) :: E
      real, dimension(ics:ice,jcs:jce,ifs:ife) :: K
      
      real, dimension(ics:ice,jcs:jce,ifs:ife) :: flux_x
      real, dimension(ics:ice,jcs:jce,ifs:ife) :: flux_y
      
      real, dimension(ics:ice,jcs:jce,ifs:ife) :: vorticity
      
      real, dimension(ids:ide,jds:jde,ifs:ife) :: flux_x_ext_p
      real, dimension(ids:ide,jds:jde,ifs:ife) :: flux_x_ext_n
      real, dimension(ids:ide,jds:jde,ifs:ife) :: flux_y_ext_p
      real, dimension(ids:ide,jds:jde,ifs:ife) :: flux_y_ext_n
      real, dimension(ids:ide,jds:jde,ifs:ife) :: E_ext_p
      real, dimension(ids:ide,jds:jde,ifs:ife) :: E_ext_n
      real, dimension(ids:ide,jds:jde,ifs:ife) :: phiG_ext_p
      real, dimension(ids:ide,jds:jde,ifs:ife) :: phiG_ext_n
      real, dimension(ids:ide,jds:jde,ifs:ife) :: u_ext_p
      real, dimension(ids:ide,jds:jde,ifs:ife) :: u_ext_n
      real, dimension(ids:ide,jds:jde,ifs:ife) :: v_ext_p
      real, dimension(ids:ide,jds:jde,ifs:ife) :: v_ext_n
      
      real, dimension(ids:ide,jds:jde,ifs:ife) :: flux_x_ext
      real, dimension(ids:ide,jds:jde,ifs:ife) :: flux_y_ext
      real, dimension(ids:ide,jds:jde,ifs:ife) :: E_ext
      real, dimension(ids:ide,jds:jde,ifs:ife) :: phiG_ext
      real, dimension(ids:ide,jds:jde,ifs:ife) :: u_ext
      real, dimension(ids:ide,jds:jde,ifs:ife) :: v_ext
      
      real, dimension(ics:ice,jcs:jce,ifs:ife) :: dfluxdx
      real, dimension(ics:ice,jcs:jce,ifs:ife) :: dfluxdy
      real, dimension(ics:ice,jcs:jce,ifs:ife) :: dEdx
      real, dimension(ics:ice,jcs:jce,ifs:ife) :: dEdy
      real, dimension(ics:ice,jcs:jce,ifs:ife) :: dvdx
      real, dimension(ics:ice,jcs:jce,ifs:ife) :: dudy
      
      integer i,j,iPatch
      
      K = 0.5 * ( stat%u * stat%uc + stat%v * stat%vc )
      E = stat%phi + mesh%phis + K
      
      flux_x = stat%phiG * stat%uc
      flux_y = stat%phiG * stat%vc
      
      do iPatch = ifs, ife
        do j = jts, jte
          do i = its, ite
            dfluxdx(i,j,iPatch) = ( flux_x(i+1,j,iPatch) - flux_x(i-1,j,iPatch) ) / (2.*dx)
            dfluxdy(i,j,iPatch) = ( flux_x(i,j+1,iPatch) - flux_x(i,j-1,iPatch) ) / (2.*dy)
            dEdx   (i,j,iPatch) = ( E     (i+1,j,iPatch) - E     (i-1,j,iPatch) ) / (2.*dx)
            dEdy   (i,j,iPatch) = ( E     (i,j+1,iPatch) - E     (i,j-1,iPatch) ) / (2.*dy)
            dvdx   (i,j,iPatch) = ( stat%v(i+1,j,iPatch) - stat%v(i-1,j,iPatch) ) / (2.*dx)
            dudy   (i,j,iPatch) = ( stat%u(i,j+1,iPatch) - stat%u(i,j-1,iPatch) ) / (2.*dy)
            
            vorticity(i,j,iPatch) = ( dvdx(i,j,iPatch) - dudy(i,j,iPatch) ) / mesh%sqrtG(i,j,iPatch) + mesh%f(i,j,iPatch)
          enddo
        enddo
      enddo
      
      tend%phiG = - ( dfluxdx + dfluxdy )
      tend%u    = - dEdx + mesh%sqrtG * vorticity * stat%vc
      tend%v    = - dEdy - mesh%sqrtG * vorticity * stat%uc
    end subroutine spatial_operator
    
    subroutine reconstruct(field_ext_p,field_ext_n,field)
      real, dimension(ids:ide,jds:jde,ifs:ife) :: field_ext_p
      real, dimension(ids:ide,jds:jde,ifs:ife) :: field_ext_n
      real, dimension(ics:ice,jcs:jce,ifs:ife) :: field
      
      integer i,j,iPatch
      
      if(trim(adjustl(reconstruct_scheme)) == 'WENO' .or. trim(adjustl(reconstruct_scheme)) == 'weno')then
        call weno(field_ext_p,field_ext_n,field)
      else
        stop 'unknown reconstruct scheme'
      endif
    
    end subroutine reconstruct
    
    subroutine weno(field_ext_p,field_ext_n,field)
      real, dimension(ids:ide,jds:jde,ifs:ife) :: field_ext_p
      real, dimension(ids:ide,jds:jde,ifs:ife) :: field_ext_n
      real, dimension(ics:ice,jcs:jce,ifs:ife) :: field
      
      real Qx(5)
      real Qy(5)
      
      integer i,j,iPatch
      integer iext
      integer jext
      
      integer dir
      
      do iPatch = ifs, ife
        do j = jts, jte
          do i = its-1, ite
            Qx(1) = field(i-2,j,iPatch)
            Qx(2) = field(i-1,j,iPatch)
            Qx(3) = field(i  ,j,iPatch)
            Qx(4) = field(i+1,j,iPatch)
            Qx(5) = field(i+2,j,iPatch)
            
            iext = 2 * i + 1
            dir  = 1
            call WENO_limiter(field_ext_p(iext,j,iPatch),Qx,dir)
          enddo
          do i = its, ite+1
            Qx(1) = field(i-2,j,iPatch)
            Qx(2) = field(i-1,j,iPatch)
            Qx(3) = field(i  ,j,iPatch)
            Qx(4) = field(i+1,j,iPatch)
            Qx(5) = field(i+2,j,iPatch)
            
            iext = 2 * i - 1
            dir  = -1
            call WENO_limiter(field_ext_n(iext,j,iPatch),Qx,dir)
          enddo
        enddo
        
        do i = its, ite
          do j = jts-1, jte
            Qy(1) = field(i,j-2,iPatch)
            Qy(2) = field(i,j-1,iPatch)
            Qy(3) = field(i,j  ,iPatch)
            Qy(4) = field(i,j+1,iPatch)
            Qy(5) = field(i,j+2,iPatch)
            
            jext = 2 * j + 1
            dir  = 1
            call WENO_limiter(field_ext_p(i,jext,iPatch),Qy,dir)
          enddo
          do j = jts, jte+1
            Qy(1) = field(i,j-2,iPatch)
            Qy(2) = field(i,j-1,iPatch)
            Qy(3) = field(i,j  ,iPatch)
            Qy(4) = field(i,j+1,iPatch)
            Qy(5) = field(i,j+2,iPatch)
            
            jext = 2 * j - 1
            dir  = -1
            call WENO_limiter(field_ext_n(i,jext,iPatch),Qy,dir)
          enddo
        enddo
      enddo

    end subroutine weno
    
    ! 1D WENO slope limiter, according to Sun,2015
    ! "A Slope Constrained 4th OrderMulti-Moment Finite Volume Method with WENO Limiter"
    ! and Jiang and Shu, 1996
    subroutine WENO_limiter(Qrec,Q,dir)
      real              , intent(out) :: Qrec
      real, dimension(5), intent(in ) :: Q
      integer           , intent(in ) :: dir
      
      integer, parameter :: nStencil = 3
      real   , parameter :: weno_coef(3)  = [0.1, 0.6, 0.3]
      real   , parameter :: eps           = 1.E-2
      
      real Qim(nStencil-1)
      real Qip(nStencil-1)
      real Qi
      
      real, dimension(nStencil) :: stencil
      real, dimension(nStencil) :: coefA
      real, dimension(nStencil) :: coefB
      real, dimension(nStencil) :: alpha
      real, dimension(nStencil) :: beta
      real, dimension(nStencil) :: omega
      
      real tau40
      real tau41
      real tau5
      
      integer iStencil
      
      Qim(2) = Q(1)
      Qim(1) = Q(2)
      Qi     = Q(3)
      Qip(1) = Q(4)
      Qip(2) = Q(5)
      
      if(dir>0)then
        stencil (1) =  Qim(2)/3. - 7./6. * Qim(1) + 11./6. * Qi     
        stencil (2) = -Qim(1)/6. + 5./6. * Qi     +  1./3. * Qip(1) 
        stencil (3) =  Qi    /3. + 5./6. * Qip(1) -  1./6. * Qip(2)
        
        coefA(1) = Qim(2) - 2. * Qim(1) + Qi
        coefA(2) = Qim(1) - 2. * Qi     + Qip(1)
        coefA(3) = Qi     - 2. * Qip(1) + Qip(2)
        
        coefB(1) =      Qim(2) - 4. * Qim(1) + 3. * Qi
        coefB(2) =      Qim(1) -      Qip(1)
        coefB(3) = 3. * Qi     - 4. * Qip(1) +      Qip(2)
      elseif(dir<0)then
        stencil (1) =  Qip(2)/3. - 7./6. * Qip(1) + 11./6. * Qi     
        stencil (2) = -Qip(1)/6. + 5./6. * Qi     +  1./3. * Qim(1) 
        stencil (3) =  Qi    /3. + 5./6. * Qim(1) -  1./6. * Qim(2)
        
        coefA(1) = Qip(2) - 2. * Qip(1) + Qi
        coefA(2) = Qip(1) - 2. * Qi     + Qim(1)
        coefA(3) = Qi     - 2. * Qim(1) + Qim(2)
        
        coefB(1) =      Qip(2) - 4. * Qip(1) + 3. * Qi
        coefB(2) =      Qip(1) -      Qim(1)
        coefB(3) = 3. * Qi     - 4. * Qim(1) +      Qim(2)
      endif
      
      beta = coefA**2 * 13. / 12. + coefB**2 * 0.25
      
      ! WENO-Z
      tau40 = abs( beta(1) - beta(2) )
      tau41 = abs( beta(2) - beta(3) )
      tau5  = abs( beta(3) - beta(1) )
      
      ! xi
      if( tau40<=minval(beta) .and. tau41>minval(beta) )then
        omega = [1./3.,2./3.,0.]
      elseif( tau40>minval(beta) .and. tau41<minval(beta) )then
        omega = [0.,2./3.,1./3.]
      else
        alpha = weno_coef * ( 1. + tau5 / ( eps + beta ) )
        omega = alpha / sum(alpha)
      endif
      
      !omega = weno_coef
      
      Qrec = dot_product( stencil, omega )
      
    end subroutine WENO_limiter
    
    !subroutine unify_bdy_stat(stat)
    !  type(stat_field), intent(inout) :: stat
    !  
    !  call unify_bdy_field(stat%phi            (ids:ide,jds:jde,ifs:ife))
    !  call unify_bdy_field(stat%zonal_wind     (ids:ide,jds:jde,ifs:ife))
    !  call unify_bdy_field(stat%meridional_wind(ids:ide,jds:jde,ifs:ife))
    !  
    !end subroutine unify_bdy_stat
    
    subroutine unify_bdy_field(field)
      real,intent(inout) :: field(ids:ide,jds:jde,ifs:ife)
      
      real field_raw(ids:ide,jds:jde,ifs:ife)
      
      field_raw = field
      
      ! Bdy lines
      ! edge(1,4)
      field(ids,jds:jde,1) = 0.5 * (field_raw(ide,jds:jde,4) + field_raw(ids,jds:jde,1)) ! left bdy of patch 1
      field(ide,jds:jde,4) = field(ids,jds:jde,1)
      
      ! edge(1,2)
      field(ids,jds:jde,2) = 0.5 * (field_raw(ide,jds:jde,1) + field_raw(ids,jds:jde,2)) ! left bdy of patch 2
      field(ide,jds:jde,1) = field(ids,jds:jde,2)
      
      ! edge(2,3)
      field(ids,jds:jde,3) = 0.5 * (field_raw(ide,jds:jde,2) + field_raw(ids,jds:jde,3)) ! left bdy of patch 3
      field(ide,jds:jde,2) = field(ids,jds:jde,3)
      
      ! edge(3,4)
      field(ids,jds:jde,4) = 0.5 * (field_raw(ide,jds:jde,3) + field_raw(ids,jds:jde,4)) ! left bdy of patch 4
      field(ide,jds:jde,3) = field(ids,jds:jde,4)
      
      ! edge(1,5)
      field(ids:ide,jde,1) = 0.5 * (field_raw(ids:ide,jds,5) + field_raw(ids:ide,jde,1)) ! top bdy of patch 1
      field(ids:ide,jds,5) = field(ids:ide,jde,1)
      
      ! edge(1,6)
      field(ids:ide,jds,1) = 0.5 * (field_raw(ids:ide,jds,1) + field_raw(ids:ide,jde,6)) ! bottom bdy of patch 1
      field(ids:ide,jde,6) = field(ids:ide,jds,1)
      
      ! edge(2,5)
      field(ids:ide,jde,2) = 0.5 * (field_raw(ide,jds:jde,5) + field_raw(ids:ide,jde,2)) ! top bdy of patch 2
      field(ide,jds:jde,5) = field(ids:ide,jde,2)
      
      ! edge(2,6)
      field(ids:ide,jds   ,2) = 0.5 * (field_raw(ids:ide,jds,2) + field_raw(ide,jde:jds:-1,6)) ! bottom bdy of patch 2
      field(ide,jde:jds:-1,6) = field(ids:ide,jds,2)
      
      ! edge(3,5)
      field(ids:ide,jde   ,3) = 0.5 * (field_raw(ids:ide,jde,3) + field_raw(ide:ids:-1,jde,5)) ! top bdy of patch 3
      field(ide:ids:-1,jde,5) = field(ids:ide,jde,3)
      
      ! edge(3,6)
      field(ids:ide,jds   ,3) = 0.5 * (field_raw(ids:ide,jds,3) + field_raw(ide:ids:-1,jds,6)) ! bottom bdy of patch 3
      field(ide:ids:-1,jds,6) = field(ids:ide,jds,3)
      
      ! edge(4,5)
      field(ids:ide,jde   ,4) = 0.5 * (field_raw(ids:ide,jde,4) + field_raw(ids,jde:jds:-1,5)) ! top bdy of patch 4
      field(ids,jde:jds:-1,5) = field(ids:ide,jde,4)
      
      ! edge(4,6)
      field(ids:ide,jds,4) = 0.5 * (field_raw(ids:ide,jds,4) + field_raw(ids,jds:jde,6)) ! bottom bdy of patch 4
      field(ids,jds:jde,6) = field(ids:ide,jds,4)
      
    end subroutine unify_bdy_field
    
    ! convert vector from patch to sphere
    subroutine convert_wind_P2SP(stat)
      type(stat_field), target, intent(inout) :: stat
      
      integer i,j,iPatch
      
      !$OMP PARALLEL DO PRIVATE(i,j)
      do iPatch = ifs, ife
        do j = jts, jte
          do i = its, ite
            call covProjPlane2Sphere(stat%zonal_wind(i,j,iPatch), stat%meridional_wind(i,j,iPatch), stat%u(i,j,iPatch), stat%v(i,j,iPatch), mesh%matrixA(:,:,i,j,iPatch), mesh%matrixIG(:,:,i,j,iPatch))
          enddo
        enddo
      enddo
      !$OMP END PARALLEL DO
    end subroutine convert_wind_P2SP
    
    subroutine convert_wind_SP2P(stat)
      type(stat_field), intent(inout) :: stat
      
      integer i,j,iPatch
      
      !$OMP PARALLEL DO PRIVATE(i,j)
      do iPatch = ifs, ife
        do j = jcs, jce
          do i = ics, ice
            call covProjSphere2Plane(stat%u(i,j,iPatch), stat%v(i,j,iPatch), stat%zonal_wind(i,j,iPatch), stat%meridional_wind(i,j,iPatch), mesh%matrixIA(:,:,i,j,iPatch), mesh%matrixG(:,:,i,j,iPatch))
          enddo
        enddo
      enddo
      !$OMP END PARALLEL DO
    end subroutine convert_wind_SP2P
    
    ! convert vector from patch to sphere near the boundary of each patch, for interpolating ghost zone
    subroutine convert_ghost_wind_P2SP(stat)
      type(stat_field), target, intent(inout) :: stat
      
      integer i,j,iPatch
      
      !$OMP PARALLEL DO PRIVATE(i,j)
      do iPatch = ifs, ife
        ! left
        do j = jts, jte
          do i = its+1, its+xhalo
            call covProjPlane2Sphere(stat%zonal_wind(i,j,iPatch), stat%meridional_wind(i,j,iPatch), stat%u(i,j,iPatch), stat%v(i,j,iPatch), mesh%matrixA(:,:,i,j,iPatch), mesh%matrixIG(:,:,i,j,iPatch))
          enddo
        enddo
        ! right
        do j = jts, jte
          do i = ite-xhalo, ite
            call covProjPlane2Sphere(stat%zonal_wind(i,j,iPatch), stat%meridional_wind(i,j,iPatch), stat%u(i,j,iPatch), stat%v(i,j,iPatch), mesh%matrixA(:,:,i,j,iPatch), mesh%matrixIG(:,:,i,j,iPatch))
          enddo
        enddo
        ! top
        do j = jte-xhalo, jte
          do i = its+xhalo, ite-xhalo
            call covProjPlane2Sphere(stat%zonal_wind(i,j,iPatch), stat%meridional_wind(i,j,iPatch), stat%u(i,j,iPatch), stat%v(i,j,iPatch), mesh%matrixA(:,:,i,j,iPatch), mesh%matrixIG(:,:,i,j,iPatch))
          enddo
        enddo
        ! bottom
        do j = jts, jts+xhalo
          do i = its+xhalo, ite-xhalo
            call covProjPlane2Sphere(stat%zonal_wind(i,j,iPatch), stat%meridional_wind(i,j,iPatch), stat%u(i,j,iPatch), stat%v(i,j,iPatch), mesh%matrixA(:,:,i,j,iPatch), mesh%matrixIG(:,:,i,j,iPatch))
          enddo
        enddo
      enddo
      !$OMP END PARALLEL DO
    end subroutine convert_ghost_wind_P2SP
    
    subroutine convert_ghost_wind_SP2P(stat)
      type(stat_field), target, intent(inout) :: stat
      
      integer i,j,iPatch
      
      !$OMP PARALLEL DO PRIVATE(i,j)
      do iPatch = ifs, ife
        ! left
        do j = jts, jte
          do i = ics, its-1
            call covProjSphere2Plane(stat%u(i,j,iPatch), stat%v(i,j,iPatch), stat%zonal_wind(i,j,iPatch), stat%meridional_wind(i,j,iPatch), mesh%matrixIA(:,:,i,j,iPatch), mesh%matrixG(:,:,i,j,iPatch))
          enddo
        enddo
        ! right
        do j = jts, jte
          do i = ite+1, ice
            call covProjSphere2Plane(stat%u(i,j,iPatch), stat%v(i,j,iPatch), stat%zonal_wind(i,j,iPatch), stat%meridional_wind(i,j,iPatch), mesh%matrixIA(:,:,i,j,iPatch), mesh%matrixG(:,:,i,j,iPatch))
          enddo
        enddo
        ! top
        do j = jte+1, jce
          do i = its, ite
            call covProjSphere2Plane(stat%u(i,j,iPatch), stat%v(i,j,iPatch), stat%zonal_wind(i,j,iPatch), stat%meridional_wind(i,j,iPatch), mesh%matrixIA(:,:,i,j,iPatch), mesh%matrixG(:,:,i,j,iPatch))
          enddo
        enddo
        ! bottom
        do j = jcs, jts-1
          do i = its, ite
            call covProjSphere2Plane(stat%u(i,j,iPatch), stat%v(i,j,iPatch), stat%zonal_wind(i,j,iPatch), stat%meridional_wind(i,j,iPatch), mesh%matrixIA(:,:,i,j,iPatch), mesh%matrixG(:,:,i,j,iPatch))
          enddo
        enddo
      enddo
      !$OMP END PARALLEL DO
    end subroutine convert_ghost_wind_SP2P
    
    subroutine convert_wind_cov2contrav(stat)
      type(stat_field), target, intent(inout) :: stat
      
      integer i,j,iPatch
      
      !$OMP PARALLEL DO PRIVATE(i,j)
      do iPatch = ifs, ife
        do j = jcs, jce
          do i = ics, ice
            call cov2contrav(stat%uc(i,j,iPatch),stat%vc(i,j,iPatch),stat%u(i,j,iPatch),stat%v(i,j,iPatch),mesh%matrixIG(:,:,i,j,iPatch))
          enddo
        enddo
      enddo
      !$OMP END PARALLEL DO
    end subroutine convert_wind_cov2contrav
    
END MODULE spatial_operators_mod

