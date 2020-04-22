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
      
      real eigenvalue_x(2,3)
      real eigenvalue_y(2,3)
      real maxeigen_x
      real maxeigen_y
      
      integer i,j,iPatch
      integer ic,jc
      integer ip1,jp1
      integer im1,jm1
      
      K = 0.5 * ( stat%u * stat%uc + stat%v * stat%vc )
      E = stat%phi + mesh%phis + K
      
      flux_x = stat%phiG * stat%uc
      flux_y = stat%phiG * stat%vc
      
      call reconstruct(flux_x_ext_p, flux_x_ext_n, flux_x   )
      call reconstruct(flux_y_ext_p, flux_y_ext_n, flux_y   )
      call reconstruct(E_ext_p     , E_ext_n     , E        )
      call reconstruct(phiG_ext_p  , phiG_ext_n  , stat%phiG + mesh%phis * mesh%sqrtG)
      call reconstruct(u_ext_p     , u_ext_n     , stat%u   )
      call reconstruct(v_ext_p     , v_ext_n     , stat%v   )
      
      do iPatch = ifs, ife
        !$OMP PARALLEL DO PRIVATE(i,eigenvalue_x,maxeigen_x,ip1,im1,jc)
        do j = jts, jte
          do i = its-1, ite
            eigenvalue_x(1,1) = stat%uc(i-1,j,iPatch) - sqrt( mesh%matrixIG(1,1,i-1,j,iPatch) * ( stat%phi(i-1,j,iPatch) + mesh%phis(i-1,j,iPatch) ) )
            eigenvalue_x(2,1) = stat%uc(i-1,j,iPatch) + sqrt( mesh%matrixIG(1,1,i-1,j,iPatch) * ( stat%phi(i-1,j,iPatch) + mesh%phis(i-1,j,iPatch) ) )
            eigenvalue_x(1,2) = stat%uc(i  ,j,iPatch) - sqrt( mesh%matrixIG(1,1,i  ,j,iPatch) * ( stat%phi(i  ,j,iPatch) + mesh%phis(i  ,j,iPatch) ) )
            eigenvalue_x(2,2) = stat%uc(i  ,j,iPatch) + sqrt( mesh%matrixIG(1,1,i  ,j,iPatch) * ( stat%phi(i  ,j,iPatch) + mesh%phis(i  ,j,iPatch) ) )
            eigenvalue_x(1,3) = stat%uc(i+1,j,iPatch) - sqrt( mesh%matrixIG(1,1,i+1,j,iPatch) * ( stat%phi(i+1,j,iPatch) + mesh%phis(i+1,j,iPatch) ) )
            eigenvalue_x(2,3) = stat%uc(i+1,j,iPatch) + sqrt( mesh%matrixIG(1,1,i+1,j,iPatch) * ( stat%phi(i+1,j,iPatch) + mesh%phis(i+1,j,iPatch) ) )
            
            maxeigen_x = maxval(abs(eigenvalue_x))
            
            ip1 = 2 * i + 1
            im1 = 2 * i - 1
            jc  = 2 * j
            
            flux_x_ext(ip1,jc,iPatch) = 0.5 * ( flux_x_ext_p(ip1,jc,iPatch) + flux_x_ext_n(ip1,jc,iPatch) - maxeigen_x                   * ( phiG_ext_n(ip1,jc,iPatch) - phiG_ext_p(ip1,jc,iPatch) ) )
            E_ext     (ip1,jc,iPatch) = 0.5 * ( E_ext_p     (ip1,jc,iPatch) + E_ext_n     (ip1,jc,iPatch) - maxeigen_x                   * ( u_ext_n   (ip1,jc,iPatch) - u_ext_p   (ip1,jc,iPatch) ) )
            v_ext     (ip1,jc,iPatch) = 0.5 * ( v_ext_p     (ip1,jc,iPatch) + v_ext_n     (ip1,jc,iPatch) - sign(1.,stat%uc(i,j,iPatch)) * ( v_ext_n   (ip1,jc,iPatch) - v_ext_p   (ip1,jc,iPatch) ) )
          enddo
        enddo
        !$OMP END PARALLEL DO
        
        !$OMP PARALLEL DO PRIVATE(i,eigenvalue_y,maxeigen_y,jp1,jm1,ic)
        do j = jts-1, jte
          do i = its, ite
            eigenvalue_y(1,1) = stat%vc(i,j-1,iPatch) - sqrt( mesh%matrixIG(2,2,i,j-1,iPatch) * ( stat%phi(i,j-1,iPatch) + mesh%phis(i,j-1,iPatch) ) )
            eigenvalue_y(2,1) = stat%vc(i,j-1,iPatch) + sqrt( mesh%matrixIG(2,2,i,j-1,iPatch) * ( stat%phi(i,j-1,iPatch) + mesh%phis(i,j-1,iPatch) ) )
            eigenvalue_y(1,2) = stat%vc(i,j  ,iPatch) - sqrt( mesh%matrixIG(2,2,i,j  ,iPatch) * ( stat%phi(i,j  ,iPatch) + mesh%phis(i,j  ,iPatch) ) )
            eigenvalue_y(2,2) = stat%vc(i,j  ,iPatch) + sqrt( mesh%matrixIG(2,2,i,j  ,iPatch) * ( stat%phi(i,j  ,iPatch) + mesh%phis(i,j  ,iPatch) ) )
            eigenvalue_y(1,3) = stat%vc(i,j+1,iPatch) - sqrt( mesh%matrixIG(2,2,i,j+1,iPatch) * ( stat%phi(i,j+1,iPatch) + mesh%phis(i,j+1,iPatch) ) )
            eigenvalue_y(2,3) = stat%vc(i,j+1,iPatch) + sqrt( mesh%matrixIG(2,2,i,j+1,iPatch) * ( stat%phi(i,j+1,iPatch) + mesh%phis(i,j+1,iPatch) ) )
            
            maxeigen_y = maxval(abs(eigenvalue_y))
            
            ic  = 2 * i
            jp1 = 2 * j + 1
            jm1 = 2 * j - 1
            
            flux_y_ext(ic,jp1,iPatch) = 0.5 * ( flux_y_ext_p(ic,jp1,iPatch) + flux_y_ext_n(ic,jp1,iPatch) - maxeigen_y                   * ( phiG_ext_n(ic,jp1,iPatch) - phiG_ext_p(ic,jp1,iPatch) ) )
            E_ext     (ic,jp1,iPatch) = 0.5 * ( E_ext_p     (ic,jp1,iPatch) + E_ext_n     (ic,jp1,iPatch) - maxeigen_y                   * ( v_ext_n   (ic,jp1,iPatch) - v_ext_p   (ic,jp1,iPatch) ) )
            u_ext     (ic,jp1,iPatch) = 0.5 * ( u_ext_p     (ic,jp1,iPatch) + u_ext_n     (ic,jp1,iPatch) - sign(1.,stat%vc(i,j,iPatch)) * ( u_ext_n   (ic,jp1,iPatch) - u_ext_p   (ic,jp1,iPatch) ) )
          enddo
        enddo
        !$OMP END PARALLEL DO
        
        !$OMP PARALLEL DO PRIVATE(i,ic,jc,ip1,im1,jp1,jm1)
        do j = jts, jte
          do i = its, ite
            ic  = 2 * i
            jc  = 2 * j
            ip1 = 2 * i + 1
            im1 = 2 * i - 1
            jp1 = 2 * j + 1
            jm1 = 2 * j - 1
            
            dfluxdx  (i,j,iPatch) = ( flux_x_ext(ip1,jc,iPatch) - flux_x_ext(im1,jc,iPatch) ) / dx
            dEdx     (i,j,iPatch) = ( E_ext     (ip1,jc,iPatch) - E_ext     (im1,jc,iPatch) ) / dx
            dvdx     (i,j,iPatch) = ( v_ext     (ip1,jc,iPatch) - v_ext     (im1,jc,iPatch) ) / dx
            dfluxdy  (i,j,iPatch) = ( flux_y_ext(ic,jp1,iPatch) - flux_y_ext(ic,jm1,iPatch) ) / dy
            dEdy     (i,j,iPatch) = ( E_ext     (ic,jp1,iPatch) - E_ext     (ic,jm1,iPatch) ) / dy
            dudy     (i,j,iPatch) = ( u_ext     (ic,jp1,iPatch) - u_ext     (ic,jm1,iPatch) ) / dy
            
            vorticity(i,j,iPatch) = dvdx(i,j,iPatch) - dudy(i,j,iPatch) + mesh%sqrtG(i,j,iPatch) * mesh%f(i,j,iPatch)
          enddo
        enddo
        !$OMP END PARALLEL DO
      enddo
      
      tend%phiG = - ( dfluxdx + dfluxdy )
      tend%u    = - dEdx + vorticity * stat%vc
      tend%v    = - dEdy - vorticity * stat%uc
      
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
      integer ic,jc
      integer iext
      integer jext
      
      integer dir
      
      do iPatch = ifs, ife
        !$OMP PARALLEL DO PRIVATE(i,Qx,jc,iext,dir)
        do j = jts, jte
          dir  = 1
          do i = its-1, ite
            Qx(1) = field(i-2,j,iPatch)
            Qx(2) = field(i-1,j,iPatch)
            Qx(3) = field(i  ,j,iPatch)
            Qx(4) = field(i+1,j,iPatch)
            Qx(5) = field(i+2,j,iPatch)
            
            jc   = 2 * j
            iext = 2 * i + 1
            call WENO_limiter(field_ext_p(iext,jc,iPatch),Qx,dir)
          enddo
          
          dir  = -1
          do i = its, ite+1
            Qx(1) = field(i-2,j,iPatch)
            Qx(2) = field(i-1,j,iPatch)
            Qx(3) = field(i  ,j,iPatch)
            Qx(4) = field(i+1,j,iPatch)
            Qx(5) = field(i+2,j,iPatch)
            
            jc   = 2 * j
            iext = 2 * i - 1
            call WENO_limiter(field_ext_n(iext,jc,iPatch),Qx,dir)
          enddo
        enddo
        !$OMP END PARALLEL DO
        
        !$OMP PARALLEL DO PRIVATE(j,Qy,ic,jext,dir)
        do i = its, ite
          dir  = 1
          do j = jts-1, jte
            Qy(1) = field(i,j-2,iPatch)
            Qy(2) = field(i,j-1,iPatch)
            Qy(3) = field(i,j  ,iPatch)
            Qy(4) = field(i,j+1,iPatch)
            Qy(5) = field(i,j+2,iPatch)
            
            ic   = 2 * i
            jext = 2 * j + 1
            call WENO_limiter(field_ext_p(ic,jext,iPatch),Qy,dir)
          enddo
          
          dir  = -1
          do j = jts, jte+1
            Qy(1) = field(i,j-2,iPatch)
            Qy(2) = field(i,j-1,iPatch)
            Qy(3) = field(i,j  ,iPatch)
            Qy(4) = field(i,j+1,iPatch)
            Qy(5) = field(i,j+2,iPatch)
            
            ic   = 2 * i
            jext = 2 * j - 1
            call WENO_limiter(field_ext_n(ic,jext,iPatch),Qy,dir)
          enddo
        enddo
        !$OMP END PARALLEL DO
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

