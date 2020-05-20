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
      
      real, dimension(ics:ice,jcs:jce,ifs:ife) :: phit
      real, dimension(ics:ice,jcs:jce,ifs:ife) :: phitG
      
      real, dimension(ics:ice,jcs:jce,ifs:ife) :: dfluxdx
      real, dimension(ics:ice,jcs:jce,ifs:ife) :: dfluxdy
      real, dimension(ics:ice,jcs:jce,ifs:ife) :: dEdx
      real, dimension(ics:ice,jcs:jce,ifs:ife) :: dEdy
      real, dimension(ics:ice,jcs:jce,ifs:ife) :: dvdx
      real, dimension(ics:ice,jcs:jce,ifs:ife) :: dudy
      
      real, dimension(ids:ide,jds:jde,ifs:ife) :: flux_x_ext_p
      real, dimension(ids:ide,jds:jde,ifs:ife) :: flux_x_ext_n
      real, dimension(ids:ide,jds:jde,ifs:ife) :: flux_y_ext_p
      real, dimension(ids:ide,jds:jde,ifs:ife) :: flux_y_ext_n
      real, dimension(ids:ide,jds:jde,ifs:ife) :: E_ext_p
      real, dimension(ids:ide,jds:jde,ifs:ife) :: E_ext_n
      real, dimension(ids:ide,jds:jde,ifs:ife) :: phiG_ext_p
      real, dimension(ids:ide,jds:jde,ifs:ife) :: phiG_ext_n
      real, dimension(ids:ide,jds:jde,ifs:ife) :: phitG_ext_p
      real, dimension(ids:ide,jds:jde,ifs:ife) :: phitG_ext_n
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
      
      real uBdy    ! variables on cell boundary
      real vBdy    ! variables on cell boundary
      real phitBdy ! variables on cell boundary
      
      real eigenvalue_x(2)
      real eigenvalue_y(2)
      real maxeigen_x
      real maxeigen_y
      
      integer i,j,iPatch
      integer ic,jc
      integer ip1,jp1
      integer im1,jm1
      
      phit  = stat%phi + mesh%phis
      phitG = phit * mesh%sqrtG
      
      K = 0.5 * ( stat%u * stat%uc + stat%v * stat%vc )
      E = phit + K
      
      flux_x = stat%phiG * stat%uc
      flux_y = stat%phiG * stat%vc
      
      call weno(flux_x_ext_p, flux_x   , 1, 1)
      call weno(flux_x_ext_n, flux_x   ,-1, 1)
      call weno(E_ext_p     , E        , 1, 1)
      call weno(E_ext_n     , E        ,-1, 1)
      !call weno(phiG_ext_p  , stat%phiG, 1, 1)
      !call weno(phiG_ext_n  , stat%phiG,-1, 1)
      call weno(phitG_ext_p , phitG    , 1, 1)
      call weno(phitG_ext_n , phitG    ,-1, 1)
      call weno(v_ext_p     , stat%v   , 1, 1)
      call weno(v_ext_n     , stat%v   ,-1, 1)
      call weno(flux_y_ext_p, flux_y   , 1, 2)
      call weno(flux_y_ext_n, flux_y   ,-1, 2)
      call weno(E_ext_p     , E        , 1, 2)
      call weno(E_ext_n     , E        ,-1, 2)
      !call weno(phiG_ext_p  , stat%phiG, 1, 2)
      !call weno(phiG_ext_n  , stat%phiG,-1, 2)
      call weno(phitG_ext_p , phitG    , 1, 2)
      call weno(phitG_ext_n , phitG    ,-1, 2)
      call weno(u_ext_p     , stat%u   , 1, 2)
      call weno(u_ext_n     , stat%u   ,-1, 2)
      
      do iPatch = ifs, ife
        !$OMP PARALLEL DO PRIVATE(i,ip1,im1,jc,uBdy,phitBdy,eigenvalue_x,maxeigen_x)
        do j = jts, jte
          do i = its-1, ite
            ip1 = 2 * i + 1
            im1 = 2 * i - 1
            jc  = 2 * j
            
            uBdy    = 0.5 * ( stat%uc(i+1,j,iPatch) + stat%uc(i,j,iPatch) )
            phitBdy = 0.5 * ( phit   (i+1,j,iPatch) + phit   (i,j,iPatch) )
            eigenvalue_x(1) = uBdy - sqrt( mesh%matrixIG_ext(1,1,ip1,j,iPatch) * phitBdy )
            eigenvalue_x(2) = uBdy + sqrt( mesh%matrixIG_ext(1,1,ip1,j,iPatch) * phitBdy )
            
            maxeigen_x = maxval(abs(eigenvalue_x))
            
            flux_x_ext(ip1,jc,iPatch) = 0.5 * ( flux_x_ext_p(ip1,jc,iPatch) + flux_x_ext_n(ip1,jc,iPatch) - maxeigen_x    * ( phitG_ext_n (ip1,jc,iPatch) - phitG_ext_p (ip1,jc,iPatch) ) )
            E_ext     (ip1,jc,iPatch) = 0.5 * ( E_ext_p     (ip1,jc,iPatch) + E_ext_n     (ip1,jc,iPatch) - maxeigen_x    * ( u_ext_n     (ip1,jc,iPatch) - u_ext_p     (ip1,jc,iPatch) ) )
            v_ext     (ip1,jc,iPatch) = 0.5 * ( v_ext_p     (ip1,jc,iPatch) + v_ext_n     (ip1,jc,iPatch) - sign(1.,uBdy) * ( v_ext_n     (ip1,jc,iPatch) - v_ext_p     (ip1,jc,iPatch) ) )
          enddo
        enddo
        !$OMP END PARALLEL DO
      enddo
      
      do iPatch = ifs, ife
        !$OMP PARALLEL DO PRIVATE(i,jp1,jm1,ic,vBdy,phitBdy,eigenvalue_y,maxeigen_y)
        do j = jts-1, jte
          do i = its, ite
            ic  = 2 * i
            jp1 = 2 * j + 1
            jm1 = 2 * j - 1
            
            vBdy    = 0.5 * ( stat%vc(i,j+1,iPatch) + stat%vc(i,j,iPatch) )
            phitBdy = 0.5 * ( phit   (i,j+1,iPatch) + phit   (i,j,iPatch) )
            eigenvalue_y(1) = vBdy - sqrt( mesh%matrixIG_ext(2,2,i,jp1,iPatch) * phitBdy )
            eigenvalue_y(2) = vBdy + sqrt( mesh%matrixIG_ext(2,2,i,jp1,iPatch) * phitBdy )
            
            maxeigen_y = maxval(abs(eigenvalue_y))
            
            flux_y_ext(ic,jp1,iPatch) = 0.5 * ( flux_y_ext_p(ic,jp1,iPatch) + flux_y_ext_n(ic,jp1,iPatch) - maxeigen_y    * ( phitG_ext_n (ic,jp1,iPatch) - phitG_ext_p (ic,jp1,iPatch) ) )
            E_ext     (ic,jp1,iPatch) = 0.5 * ( E_ext_p     (ic,jp1,iPatch) + E_ext_n     (ic,jp1,iPatch) - maxeigen_y    * ( v_ext_n     (ic,jp1,iPatch) - v_ext_p     (ic,jp1,iPatch) ) )
            u_ext     (ic,jp1,iPatch) = 0.5 * ( u_ext_p     (ic,jp1,iPatch) + u_ext_n     (ic,jp1,iPatch) - sign(1.,vBdy) * ( u_ext_n     (ic,jp1,iPatch) - u_ext_p     (ic,jp1,iPatch) ) )
          enddo
        enddo
        !$OMP END PARALLEL DO
      enddo
        
      do iPatch = ifs, ife
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
    
    subroutine weno(field_ext,field,flux_dir,axis_dir)
      real   , dimension(ids:ide,jds:jde,ifs:ife), intent(out) :: field_ext
      real   , dimension(ics:ice,jcs:jce,ifs:ife), intent(in ) :: field
      integer,                                     intent(in ) :: flux_dir ! 1 for positive flux, -1 for negative flux
      integer,                                     intent(in ) :: axis_dir ! 1 for x axis, 2 for y axis
      
      real Qx(5)
      real Qy(5)
      
      integer i,j,iPatch
      integer ic,jc
      integer iext
      integer jext
      
      if(axis_dir==1)then
        if(flux_dir>0)then
          do iPatch = ifs, ife
            !$OMP PARALLEL DO PRIVATE(i,Qx,jc,iext)
            do j = jts, jte
              do i = its-1, ite
                Qx(1) = field(i-2,j,iPatch)
                Qx(2) = field(i-1,j,iPatch)
                Qx(3) = field(i  ,j,iPatch)
                Qx(4) = field(i+1,j,iPatch)
                Qx(5) = field(i+2,j,iPatch)
                
                jc   = 2 * j
                iext = 2 * i + 1
                call WENO_limiter(field_ext(iext,jc,iPatch),Qx,flux_dir)
              enddo
            enddo
            !$OMP END PARALLEL DO
          enddo
        elseif(flux_dir<0)then
          do iPatch = ifs, ife
            !$OMP PARALLEL DO PRIVATE(i,Qx,jc,iext)
            do j = jts, jte
              do i = its, ite+1
                Qx(1) = field(i-2,j,iPatch)
                Qx(2) = field(i-1,j,iPatch)
                Qx(3) = field(i  ,j,iPatch)
                Qx(4) = field(i+1,j,iPatch)
                Qx(5) = field(i+2,j,iPatch)
                
                jc   = 2 * j
                iext = 2 * i - 1
                call WENO_limiter(field_ext(iext,jc,iPatch),Qx,flux_dir)
              enddo
            enddo
            !$OMP END PARALLEL DO
          enddo
        endif
      endif
        
      if(axis_dir==2)then
        if(flux_dir>0)then
          do iPatch = ifs, ife
            !$OMP PARALLEL DO PRIVATE(j,Qy,ic,jext)
            do i = its, ite
              do j = jts-1, jte
                Qy(1) = field(i,j-2,iPatch)
                Qy(2) = field(i,j-1,iPatch)
                Qy(3) = field(i,j  ,iPatch)
                Qy(4) = field(i,j+1,iPatch)
                Qy(5) = field(i,j+2,iPatch)
                
                ic   = 2 * i
                jext = 2 * j + 1
                call WENO_limiter(field_ext(ic,jext,iPatch),Qy,flux_dir)
              enddo
            enddo
            !$OMP END PARALLEL DO
          enddo
        elseif(flux_dir<0)then
          do iPatch = ifs, ife
            !$OMP PARALLEL DO PRIVATE(j,Qy,ic,jext)
            do i = its, ite
              do j = jts, jte+1
                Qy(1) = field(i,j-2,iPatch)
                Qy(2) = field(i,j-1,iPatch)
                Qy(3) = field(i,j  ,iPatch)
                Qy(4) = field(i,j+1,iPatch)
                Qy(5) = field(i,j+2,iPatch)
                
                ic   = 2 * i
                jext = 2 * j - 1
                call WENO_limiter(field_ext(ic,jext,iPatch),Qy,flux_dir)
              enddo
            enddo
            !$OMP END PARALLEL DO
          enddo
        endif
      endif
        
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
    
    ! convert vector from patch to sphere
    subroutine convert_wind_P2SP(stat)
      type(stat_field), target, intent(inout) :: stat
      
      integer i,j,iPatch
      
      do iPatch = ifs, ife
        !$OMP PARALLEL DO PRIVATE(i)
        do j = jts, jte
          do i = its, ite
            call covProjPlane2Sphere(stat%zonal_wind(i,j,iPatch), stat%meridional_wind(i,j,iPatch), stat%u(i,j,iPatch), stat%v(i,j,iPatch), mesh%matrixA(:,:,i,j,iPatch), mesh%matrixIG(:,:,i,j,iPatch))
          enddo
        enddo
        !$OMP END PARALLEL DO
      enddo
    end subroutine convert_wind_P2SP
    
    subroutine convert_wind_SP2P(stat)
      type(stat_field), intent(inout) :: stat
      
      integer i,j,iPatch
      
      do iPatch = ifs, ife
        !$OMP PARALLEL DO PRIVATE(i)
        do j = jcs, jce
          do i = ics, ice
            call covProjSphere2Plane(stat%u(i,j,iPatch), stat%v(i,j,iPatch), stat%zonal_wind(i,j,iPatch), stat%meridional_wind(i,j,iPatch), mesh%matrixIA(:,:,i,j,iPatch), mesh%matrixG(:,:,i,j,iPatch))
          enddo
        enddo
        !$OMP END PARALLEL DO
      enddo
    end subroutine convert_wind_SP2P
    
    subroutine convert_wind_cov2contrav(stat)
      type(stat_field), target, intent(inout) :: stat
      
      integer i,j,iPatch
      
      do iPatch = ifs, ife
        !$OMP PARALLEL DO PRIVATE(i)
        do j = jcs, jce
          do i = ics, ice
            call cov2contrav(stat%uc(i,j,iPatch),stat%vc(i,j,iPatch),stat%u(i,j,iPatch),stat%v(i,j,iPatch),mesh%matrixIG(:,:,i,j,iPatch))
          enddo
        enddo
        !$OMP END PARALLEL DO
      enddo
    end subroutine convert_wind_cov2contrav
    
END MODULE spatial_operators_mod

