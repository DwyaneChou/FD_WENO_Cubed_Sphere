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
      
      real, dimension(ics:ice,jcs:jce,ifs:ife) :: dfluxdx
      real, dimension(ics:ice,jcs:jce,ifs:ife) :: dfluxdy
      real, dimension(ics:ice,jcs:jce,ifs:ife) :: dEdx
      real, dimension(ics:ice,jcs:jce,ifs:ife) :: dEdy
      real, dimension(ics:ice,jcs:jce,ifs:ife) :: dvdx
      real, dimension(ics:ice,jcs:jce,ifs:ife) :: dudy
      
      real, dimension(4,ids:ide) :: flux_x_bdy
      real, dimension(4,ids:ide) :: flux_y_bdy
      real phiGus
      real phiGvs
      
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
      !call weno(flux_y_ext_p, flux_y   , 1, 1) ! Just for unifying bdy
      !call weno(flux_y_ext_n, flux_y   ,-1, 1) ! Just for unifying bdy
      call weno(E_ext_p     , E        , 1, 1)
      call weno(E_ext_n     , E        ,-1, 1)
      call weno(phiG_ext_p  , stat%phiG, 1, 1)
      call weno(phiG_ext_n  , stat%phiG,-1, 1)
      call weno(phitG_ext_p , phitG    , 1, 1)
      call weno(phitG_ext_n , phitG    ,-1, 1)
      call weno(u_ext_p     , stat%u   , 1, 1)
      call weno(u_ext_n     , stat%u   ,-1, 1)
      call weno(v_ext_p     , stat%v   , 1, 1)
      call weno(v_ext_n     , stat%v   ,-1, 1)
      !call weno(flux_x_ext_p, flux_x   , 1, 2) ! Just for unifying bdy
      !call weno(flux_x_ext_n, flux_x   ,-1, 2) ! Just for unifying bdy
      call weno(flux_y_ext_p, flux_y   , 1, 2)
      call weno(flux_y_ext_n, flux_y   ,-1, 2)
      call weno(E_ext_p     , E        , 1, 2)
      call weno(E_ext_n     , E        ,-1, 2)
      call weno(phiG_ext_p  , stat%phiG, 1, 2)
      call weno(phiG_ext_n  , stat%phiG,-1, 2)
      call weno(phitG_ext_p , phitG    , 1, 2)
      call weno(phitG_ext_n , phitG    ,-1, 2)
      call weno(u_ext_p     , stat%u   , 1, 2)
      call weno(u_ext_n     , stat%u   ,-1, 2)
      call weno(v_ext_p     , stat%v   , 1, 2)
      call weno(v_ext_n     , stat%v   ,-1, 2)
      
      do iPatch = ifs, ife
        !$OMP PARALLEL DO PRIVATE(i,uBdy,phitBdy,eigenvalue_x,maxeigen_x,ip1,im1,jc)
        do j = jts, jte
          do i = its-1, ite
            ip1 = 2 * i + 1
            im1 = 2 * i - 1
            jc  = 2 * j
            
            uBdy    = 0.5 * ( stat%uc(i+1,j,iPatch) + stat%uc(i,j,iPatch) )
            phitBdy = 0.5 * ( phit   (i+1,j,iPatch) + phit   (i,j,iPatch) )
            eigenvalue_x(1) = uBdy - sqrt( mesh%matrixIG_ext(1,1,ip1,j,iPatch) * phitBdy )
            eigenvalue_x(1) = uBdy + sqrt( mesh%matrixIG_ext(1,1,ip1,j,iPatch) * phitBdy )
            
            maxeigen_x = maxval(abs(eigenvalue_x))
            
            flux_x_ext(ip1,jc,iPatch) = 0.5 * ( flux_x_ext_p(ip1,jc,iPatch) + flux_x_ext_n(ip1,jc,iPatch) - maxeigen_x    * ( phitG_ext_n (ip1,jc,iPatch) - phitG_ext_p (ip1,jc,iPatch) ) )
            !flux_y_ext(ip1,jc,iPatch) = 0.5 * ( flux_y_ext_p(ip1,jc,iPatch) + flux_y_ext_n(ip1,jc,iPatch) - sign(1.,uBdy) * ( flux_y_ext_n(ip1,jc,iPatch) - flux_y_ext_p(ip1,jc,iPatch) ) ) ! Just for unifying bdy
            E_ext     (ip1,jc,iPatch) = 0.5 * ( E_ext_p     (ip1,jc,iPatch) + E_ext_n     (ip1,jc,iPatch) - maxeigen_x    * ( u_ext_n     (ip1,jc,iPatch) - u_ext_p     (ip1,jc,iPatch) ) )
            phiG_ext  (ip1,jc,iPatch) = 0.5 * ( phiG_ext_p  (ip1,jc,iPatch) + phiG_ext_n  (ip1,jc,iPatch) - sign(1.,uBdy) * ( phiG_ext_n  (ip1,jc,iPatch) - phiG_ext_p  (ip1,jc,iPatch) ) )
            u_ext     (ip1,jc,iPatch) = 0.5 * ( u_ext_p     (ip1,jc,iPatch) + u_ext_n     (ip1,jc,iPatch) - sign(1.,uBdy) * ( u_ext_n     (ip1,jc,iPatch) - u_ext_p     (ip1,jc,iPatch) ) ) ! Just for unifying bdy
            v_ext     (ip1,jc,iPatch) = 0.5 * ( v_ext_p     (ip1,jc,iPatch) + v_ext_n     (ip1,jc,iPatch) - sign(1.,uBdy) * ( v_ext_n     (ip1,jc,iPatch) - v_ext_p     (ip1,jc,iPatch) ) )
          enddo
        enddo
        !$OMP END PARALLEL DO
      enddo
      
      do iPatch = ifs, ife
        !$OMP PARALLEL DO PRIVATE(i,vBdy,phitBdy,eigenvalue_y,maxeigen_y,jp1,jm1,ic)
        do j = jts-1, jte
          do i = its, ite
            ic  = 2 * i
            jp1 = 2 * j + 1
            jm1 = 2 * j - 1
            
            vBdy    = 0.5 * ( stat%vc(i,j+1,iPatch) + stat%vc(i,j,iPatch) )
            phitBdy = 0.5 * ( phit   (i,j+1,iPatch) + phit   (i,j,iPatch) )
            eigenvalue_y(1) = vBdy - sqrt( mesh%matrixIG_ext(2,2,i,jp1,iPatch) * phitBdy )
            eigenvalue_y(1) = vBdy + sqrt( mesh%matrixIG_ext(2,2,i,jp1,iPatch) * phitBdy )
            
            maxeigen_y = maxval(abs(eigenvalue_y))
            
            !flux_x_ext(ic,jp1,iPatch) = 0.5 * ( flux_x_ext_p(ic,jp1,iPatch) + flux_x_ext_n(ic,jp1,iPatch) - sign(1.,vBdy) * ( flux_x_ext_n(ic,jp1,iPatch) - flux_x_ext_p(ic,jp1,iPatch) ) ) ! Just for unifying bdy
            flux_y_ext(ic,jp1,iPatch) = 0.5 * ( flux_y_ext_p(ic,jp1,iPatch) + flux_y_ext_n(ic,jp1,iPatch) - maxeigen_y    * ( phitG_ext_n (ic,jp1,iPatch) - phitG_ext_p (ic,jp1,iPatch) ) )
            E_ext     (ic,jp1,iPatch) = 0.5 * ( E_ext_p     (ic,jp1,iPatch) + E_ext_n     (ic,jp1,iPatch) - maxeigen_y    * ( v_ext_n     (ic,jp1,iPatch) - v_ext_p     (ic,jp1,iPatch) ) )
            phiG_ext  (ic,jp1,iPatch) = 0.5 * ( phiG_ext_p  (ic,jp1,iPatch) + phiG_ext_n  (ic,jp1,iPatch) - sign(1.,vBdy) * ( phiG_ext_n  (ic,jp1,iPatch) - phiG_ext_p  (ic,jp1,iPatch) ) )
            u_ext     (ic,jp1,iPatch) = 0.5 * ( u_ext_p     (ic,jp1,iPatch) + u_ext_n     (ic,jp1,iPatch) - sign(1.,vBdy) * ( u_ext_n     (ic,jp1,iPatch) - u_ext_p     (ic,jp1,iPatch) ) )
            v_ext     (ic,jp1,iPatch) = 0.5 * ( v_ext_p     (ic,jp1,iPatch) + v_ext_n     (ic,jp1,iPatch) - sign(1.,vBdy) * ( v_ext_n     (ic,jp1,iPatch) - v_ext_p     (ic,jp1,iPatch) ) ) ! Just for unifying bdy
          enddo
        enddo
        !$OMP END PARALLEL DO
      enddo
      
      !call unify_bdy(flux_x_ext,flux_y_ext,E_ext,phiG_ext,u_ext,v_ext)
      call unify_bdy_field_contravariant(flux_x_ext,flux_y_ext)
      !call unify_bdy_field_covariant(u_ext,v_ext)
      !call unify_bdy_field_scalar(E_ext)
        
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
    
    subroutine unify_bdy(flux_x,flux_y,E_ext,phiG_ext,u_ext,v_ext)
      real,intent(inout) :: flux_x  (ids:ide,jds:jde,ifs:ife)
      real,intent(inout) :: flux_y  (ids:ide,jds:jde,ifs:ife)
      real,intent(inout) :: E_ext   (ids:ide,jds:jde,ifs:ife)
      real,intent(inout) :: phiG_ext(ids:ide,jds:jde,ifs:ife)
      real,intent(inout) :: u_ext   (ids:ide,jds:jde,ifs:ife)
      real,intent(inout) :: v_ext   (ids:ide,jds:jde,ifs:ife)
      
      real phi  (ids:ide)
      real phi1 (ids:ide),phi2 (ids:ide)
      real phiG1(ids:ide),phiG2(ids:ide)
      
      real u  (ids:ide),v  (ids:ide)
      real u1 (ids:ide),v1 (ids:ide)
      real u2 (ids:ide),v2 (ids:ide)
      real fx (ids:ide),fy (ids:ide)
      real fx1(ids:ide),fy1(ids:ide)
      real fx2(ids:ide),fy2(ids:ide)
      real uc (ids:ide),vc (ids:ide)
      real uc1(ids:ide),vc1(ids:ide)
      real uc2(ids:ide),vc2(ids:ide)
      real us (ids:ide),vs (ids:ide)
      real us1(ids:ide),vs1(ids:ide)
      real us2(ids:ide),vs2(ids:ide)
      real E  (ids:ide)
      real E1 (ids:ide),E2 (ids:ide)
      
      integer i,j,iPatch
      integer ir,jr
      integer iPatch1, iPatch2
      
      real eig(2,ids:ide)
      real eigmax(ids:ide)
      
      ! edge(4,5)
      iPatch1 = 4
      j = jde
      do i = ids, ide
        E1 (i) = E_ext (i,j,iPatch1)
        fy1(i) = flux_y(i,j,iPatch1)
        vc1(i) = flux_y(i,j,iPatch1) / phiG_ext(i,j,iPatch1)
        
        phi1 (i) = phiG_ext(i,j,iPatch1) / mesh%sqrtG_ext(i,j,iPatch1)
        phiG1(i) = phiG_ext(i,j,iPatch1)
        u1   (i) = u_ext   (i,j,iPatch1)
        v1   (i) = v_ext   (i,j,iPatch1)
      enddo
      
      iPatch2 = 5
      i = ids
      do j = jds, jde
        E2 (j) = E_ext (i,j,iPatch2)
        fx2(j) = flux_x(i,j,iPatch2)
        uc2(j) = flux_x(i,j,iPatch2) / phiG_ext(i,j,iPatch2)
        
        phi2 (j) = phiG_ext(i,j,iPatch2) / mesh%sqrtG_ext(i,j,iPatch2)
        phiG2(j) = phiG_ext(i,j,iPatch2)
        u2   (j) = u_ext   (i,j,iPatch2)
        v2   (j) = v_ext   (i,j,iPatch2)
      enddo
      
      phi = 0.5 * ( phi1 + phi2(jde:jds:-1) )
      vc  = 0.5 * ( vc1  + uc2 (jde:jds:-1) )
      
      eig(1,:) = vc - sqrt( mesh%matrixIG_ext(2,2,:,jde,iPatch1) * phi )
      eig(2,:) = vc + sqrt( mesh%matrixIG_ext(2,2,:,jde,iPatch1) * phi )
      
      do i = ids, ide
        eigmax(i) = max(abs(eig(1,i)),abs(eig(2,i)))
      enddo
      
      ! Convert wind from patch2 to patch1
      i = ids
      do j = jds, jde
        call covProjPlane2Sphere(us2(j), vs2(j), u2(j), v2(j), mesh%matrixA_ext(:,:,i,j,iPatch2), mesh%matrixIG_ext(:,:,i,j,iPatch2))
      enddo
      
      j = jde
      do i = ids, ide
        ir = ide - i + 1
        call covProjSphere2Plane(u2(ir), v2(ir), us2(ir), vs2(ir), mesh%matrixIA_ext(:,:,i,j,iPatch1), mesh%matrixG_ext(:,:,i,j,iPatch1))
      enddo
      
      ! Calculate Riemann solver
      fx = 0.5 * ( fy1 + fx2(jde:jds:-1) - eigmax * ( phiG2(jde:jds:-1) - phiG1 ) )
      E  = 0.5 * ( E1  + E2 (jde:jds:-1) - eigmax * ( v2   (jde:jds:-1) - v1    ) )
      
      !print*,'1',E     (ids+1:ide:2 )
      !print*,'2',E1    (ids+1:ide:2 )
      !print*,'3',E2    (jde-1:jds:-2)
      !print*,'4',eigmax(ids+1:ide:2 ) * ( u2   (jde-1:jds:-2) - v1(ids+1:ide:2)    )
      !print*,'5',v1    (ids+1:ide:2 )
      !print*,'6',u2    (jde-1:jds:-2)
      !print*,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
      
      j = jde
      do i = ids, ide
        flux_y(i,j,iPatch1) = fx(i)
        E_ext (i,j,iPatch1) = E (i)
      enddo
      
      i = ids
      do j = jds, jde
        jr = jde - j + 1
        flux_x(i,j,iPatch2) = fx(jr)
        E_ext (i,j,iPatch2) = E (jr)
      enddo
      
      ! edge(1,5)
      iPatch1 = 1
      j = jde
      do i = ids, ide
        E1 (i) = E_ext (i,j,iPatch1)
        fy1(i) = flux_y(i,j,iPatch1)
        vc1(i) = flux_y(i,j,iPatch1) / phiG_ext(i,j,iPatch1)
        
        phi1 (i) = phiG_ext(i,j,iPatch1) / mesh%sqrtG_ext(i,j,iPatch1)
        phiG1(i) = phiG_ext(i,j,iPatch1)
        u1   (i) = u_ext   (i,j,iPatch1)
        v1   (i) = v_ext   (i,j,iPatch1)
      enddo
      
      iPatch2 = 5
      j = jds
      do i = ids, ide
        E2 (i) = E_ext (i,j,iPatch2)
        fy2(i) = flux_y(i,j,iPatch2)
        vc2(i) = flux_y(i,j,iPatch2) / phiG_ext(i,j,iPatch2)
        
        phi2 (i) = phiG_ext(i,j,iPatch2) / mesh%sqrtG_ext(i,j,iPatch2)
        phiG2(i) = phiG_ext(i,j,iPatch2)
        u2   (i) = u_ext   (i,j,iPatch2)
        v2   (i) = v_ext   (i,j,iPatch2)
      enddo
      
      phi = 0.5 * ( phi1 + phi2 )
      vc  = 0.5 * ( vc1  + vc2  )
      
      eig(1,:) = vc - sqrt( mesh%matrixIG_ext(2,2,:,jde,iPatch1) * phi )
      eig(2,:) = vc + sqrt( mesh%matrixIG_ext(2,2,:,jde,iPatch1) * phi )
      
      do i = ids, ide
        eigmax(i) = max(abs(eig(1,i)),abs(eig(2,i)))
      enddo
      
      ! Convert wind from patch2 to patch1
      j = jds
      do i = ids, ide
        call covProjPlane2Sphere(us2(i), vs2(i), u2(i), v2(i), mesh%matrixA_ext(:,:,i,j,iPatch2), mesh%matrixIG_ext(:,:,i,j,iPatch2))
      enddo
      
      j = jde
      do i = ids, ide
        call covProjSphere2Plane(u2(i), v2(i), us2(i), vs2(i), mesh%matrixIA_ext(:,:,i,j,iPatch1), mesh%matrixG_ext(:,:,i,j,iPatch1))
      enddo
      
      ! Calculate Riemann solver
      fx = 0.5 * ( fy1 + fy2 - eigmax * ( phiG2 - phiG1 ) )
      E  = 0.5 * ( E1  + E2  - eigmax * ( v2    - v1    ) )
      
      !print*,'1',E     (ids+1:ide:2 )
      !print*,'2',E1    (ids+1:ide:2 )
      !print*,'3',E2    (ids+1:ide:2 )
      !print*,'4',eigmax(ids+1:ide:2 ) * ( v2   (ids+1:ide:2 ) - v1(ids+1:ide:2)    )
      !print*,'5',v1    (ids+1:ide:2 )
      !print*,'6',v2    (ids+1:ide:2 )
      !print*,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
      !stop
      
      j = jde
      do i = ids, ide
        flux_y(i,j,iPatch1) = fx(i)
        E_ext (i,j,iPatch1) = E (i)
      enddo
      
      j = jds
      do i = ids, ide
        flux_y(i,j,iPatch2) = fx(i)
        E_ext (i,j,iPatch2) = E (i)
      enddo
      
    end subroutine unify_bdy
    
    subroutine unify_bdy_field_scalar(field)
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
      
    end subroutine unify_bdy_field_scalar
    
    subroutine unify_bdy_field_contravariant(field_x,field_y)
      real,intent(inout) :: field_x(ids:ide,jds:jde,ifs:ife)
      real,intent(inout) :: field_y(ids:ide,jds:jde,ifs:ife)
      
      real us (ids:ide),vs (ids:ide)
      real us1(ids:ide),vs1(ids:ide)
      real us2(ids:ide),vs2(ids:ide)
      
      integer i,j,iPatch
      integer ir,jr
      integer iPatch1, iPatch2
      
      real field_x_raw(ids:ide,jds:jde,ifs:ife)
      real field_y_raw(ids:ide,jds:jde,ifs:ife)
      
      field_x_raw = field_x
      field_y_raw = field_y
      
      ! Bdy lines
      ! edge(1,4)
      iPatch1 = 1
      i = ids
      do j = jds, jde
        call contravProjPlane2Sphere(us1(j), vs1(j), field_x_raw(i,j,iPatch1), field_y_raw(i,j,iPatch1), mesh%matrixA_ext(:,:,i,j,iPatch1))
      enddo
      
      iPatch2 = 4
      i = ide
      do j = jds, jde
        call contravProjPlane2Sphere(us2(j), vs2(j), field_x_raw(i,j,iPatch2), field_y_raw(i,j,iPatch2), mesh%matrixA_ext(:,:,i,j,iPatch2))
      enddo
      us = 0.5 * ( us1 + us2 )
      vs = 0.5 * ( vs1 + vs2 )
      
      i = ids
      do j = jds, jde
        call contravProjSphere2Plane(field_x(i,j,iPatch1), field_y(i,j,iPatch1), us(j), vs(j), mesh%matrixIA_ext(:,:,i,j,iPatch1))
      enddo
      
      i = ide
      do j = jds, jde
        call contravProjSphere2Plane(field_x(i,j,iPatch2), field_y(i,j,iPatch2), us(j), vs(j), mesh%matrixIA_ext(:,:,i,j,iPatch2))
      enddo

      ! edge(1,2)
      iPatch1 = 2
      i = ids
      do j = jds, jde
        call contravProjPlane2Sphere(us1(j), vs1(j), field_x_raw(i,j,iPatch1), field_y_raw(i,j,iPatch1), mesh%matrixA_ext(:,:,i,j,iPatch1))
      enddo
      
      iPatch2 = 1
      i = ide
      do j = jds, jde
        call contravProjPlane2Sphere(us2(j), vs2(j), field_x_raw(i,j,iPatch2), field_y_raw(i,j,iPatch2), mesh%matrixA_ext(:,:,i,j,iPatch2))
      enddo
      us = 0.5 * ( us1 + us2 )
      vs = 0.5 * ( vs1 + vs2 )
      
      i = ids
      do j = jds, jde
        call contravProjSphere2Plane(field_x(i,j,iPatch1), field_y(i,j,iPatch1), us(j), vs(j), mesh%matrixIA_ext(:,:,i,j,iPatch1))
      enddo
      
      i = ide
      do j = jds, jde
        call contravProjSphere2Plane(field_x(i,j,iPatch2), field_y(i,j,iPatch2), us(j), vs(j), mesh%matrixIA_ext(:,:,i,j,iPatch2))
      enddo
      
      ! edge(2,3)
      iPatch1 = 3
      i = ids
      do j = jds, jde
        call contravProjPlane2Sphere(us1(j), vs1(j), field_x_raw(i,j,iPatch1), field_y_raw(i,j,iPatch1), mesh%matrixA_ext(:,:,i,j,iPatch1))
      enddo
      
      iPatch2 = 2
      i = ide
      do j = jds, jde
        call contravProjPlane2Sphere(us2(j), vs2(j), field_x_raw(i,j,iPatch2), field_y_raw(i,j,iPatch2), mesh%matrixA_ext(:,:,i,j,iPatch2))
      enddo
      us = 0.5 * ( us1 + us2 )
      vs = 0.5 * ( vs1 + vs2 )
      
      i = ids
      do j = jds, jde
        call contravProjSphere2Plane(field_x(i,j,iPatch1), field_y(i,j,iPatch1), us(j), vs(j), mesh%matrixIA_ext(:,:,i,j,iPatch1))
      enddo
      
      i = ide
      do j = jds, jde
        call contravProjSphere2Plane(field_x(i,j,iPatch2), field_y(i,j,iPatch2), us(j), vs(j), mesh%matrixIA_ext(:,:,i,j,iPatch2))
      enddo
      
      ! edge(3,4)
      iPatch1 = 4
      i = ids
      do j = jds, jde
        call contravProjPlane2Sphere(us1(j), vs1(j), field_x_raw(i,j,iPatch1), field_y_raw(i,j,iPatch1), mesh%matrixA_ext(:,:,i,j,iPatch1))
      enddo
      
      iPatch2 = 3
      i = ide
      do j = jds, jde
        call contravProjPlane2Sphere(us2(j), vs2(j), field_x_raw(i,j,iPatch2), field_y_raw(i,j,iPatch2), mesh%matrixA_ext(:,:,i,j,iPatch2))
      enddo
      us = 0.5 * ( us1 + us2 )
      vs = 0.5 * ( vs1 + vs2 )
      
      i = ids
      do j = jds, jde
        call contravProjSphere2Plane(field_x(i,j,iPatch1), field_y(i,j,iPatch1), us(j), vs(j), mesh%matrixIA_ext(:,:,i,j,iPatch1))
      enddo
      
      i = ide
      do j = jds, jde
        call contravProjSphere2Plane(field_x(i,j,iPatch2), field_y(i,j,iPatch2), us(j), vs(j), mesh%matrixIA_ext(:,:,i,j,iPatch2))
      enddo
      
      ! edge(1,5)
      iPatch1 = 5
      j = jds
      do i = ids, ide
        call contravProjPlane2Sphere(us1(i), vs1(i), field_x_raw(i,j,iPatch1), field_y_raw(i,j,iPatch1), mesh%matrixA_ext(:,:,i,j,iPatch1))
      enddo
      
      iPatch2 = 1
      j = jde
      do i = ids, ide
        call contravProjPlane2Sphere(us2(i), vs2(i), field_x_raw(i,j,iPatch2), field_y_raw(i,j,iPatch2), mesh%matrixA_ext(:,:,i,j,iPatch2))
      enddo
      us = 0.5 * ( us1 + us2 )
      vs = 0.5 * ( vs1 + vs2 )
      
      j = jds
      do i = ids, ide
        call contravProjSphere2Plane(field_x(i,j,iPatch1), field_y(i,j,iPatch1), us(i), vs(i), mesh%matrixIA_ext(:,:,i,j,iPatch1))
      enddo
      
      j = jde
      do i = ids, ide
        call contravProjSphere2Plane(field_x(i,j,iPatch2), field_y(i,j,iPatch2), us(i), vs(i), mesh%matrixIA_ext(:,:,i,j,iPatch2))
      enddo
      
      ! edge(1,6)
      iPatch1 = 1
      j = jds
      do i = ids, ide
        call contravProjPlane2Sphere(us1(i), vs1(i), field_x_raw(i,j,iPatch1), field_y_raw(i,j,iPatch1), mesh%matrixA_ext(:,:,i,j,iPatch1))
      enddo
      
      iPatch2 = 6
      j = jde
      do i = ids, ide
        call contravProjPlane2Sphere(us2(i), vs2(i), field_x_raw(i,j,iPatch2), field_y_raw(i,j,iPatch2), mesh%matrixA_ext(:,:,i,j,iPatch2))
      enddo
      us = 0.5 * ( us1 + us2 )
      vs = 0.5 * ( vs1 + vs2 )
      
      j = jds
      do i = ids, ide
        call contravProjSphere2Plane(field_x(i,j,iPatch1), field_y(i,j,iPatch1), us(i), vs(i), mesh%matrixIA_ext(:,:,i,j,iPatch1))
      enddo
      
      j = jde
      do i = ids, ide
        call contravProjSphere2Plane(field_x(i,j,iPatch2), field_y(i,j,iPatch2), us(i), vs(i), mesh%matrixIA_ext(:,:,i,j,iPatch2))
      enddo
      
      ! edge(2,5)
      iPatch1 = 2
      j = jde
      do i = ids, ide
        call contravProjPlane2Sphere(us1(i), vs1(i), field_x_raw(i,j,iPatch1), field_y_raw(i,j,iPatch1), mesh%matrixA_ext(:,:,i,j,iPatch1))
      enddo
      
      iPatch2 = 5
      i = ide
      do j = jds, jde
        call contravProjPlane2Sphere(us2(j), vs2(j), field_x_raw(i,j,iPatch2), field_y_raw(i,j,iPatch2), mesh%matrixA_ext(:,:,i,j,iPatch2))
      enddo
      us = 0.5 * ( us1 + us2 )
      vs = 0.5 * ( vs1 + vs2 )
      
      j = jde
      do i = ids, ide
        call contravProjSphere2Plane(field_x(i,j,iPatch1), field_y(i,j,iPatch1), us(i), vs(i), mesh%matrixIA_ext(:,:,i,j,iPatch1))
      enddo
      
      i = ide
      do j = jds, jde
        call contravProjSphere2Plane(field_x(i,j,iPatch2), field_y(i,j,iPatch2), us(j), vs(j), mesh%matrixIA_ext(:,:,i,j,iPatch2))
      enddo
      
      ! edge(2,6)
      iPatch1 = 2
      j = jds
      do i = ids, ide
        call contravProjPlane2Sphere(us1(i), vs1(i), field_x_raw(i,j,iPatch1), field_y_raw(i,j,iPatch1), mesh%matrixA_ext(:,:,i,j,iPatch1))
      enddo
      
      iPatch2 = 6
      i = ide
      do j = jds, jde
        call contravProjPlane2Sphere(us2(j), vs2(j), field_x_raw(i,j,iPatch2), field_y_raw(i,j,iPatch2), mesh%matrixA_ext(:,:,i,j,iPatch2))
      enddo
      us = 0.5 * ( us1 + us2(jde:jds:-1) )
      vs = 0.5 * ( vs1 + vs2(jde:jds:-1) )
      
      j = jds
      do i = ids, ide
        call contravProjSphere2Plane(field_x(i,j,iPatch1), field_y(i,j,iPatch1), us(i), vs(i), mesh%matrixIA_ext(:,:,i,j,iPatch1))
      enddo
      
      i = ide
      do j = jds, jde
        jr = jde - j + 1
        call contravProjSphere2Plane(field_x(i,j,iPatch2), field_y(i,j,iPatch2), us(jr), vs(jr), mesh%matrixIA_ext(:,:,i,j,iPatch2))
      enddo
      
      ! edge(3,5)
      iPatch1 = 3
      j = jde
      do i = ids, ide
        call contravProjPlane2Sphere(us1(i), vs1(i), field_x_raw(i,j,iPatch1), field_y_raw(i,j,iPatch1), mesh%matrixA_ext(:,:,i,j,iPatch1))
      enddo
      
      iPatch2 = 5
      j = jde
      do i = ids, ide
        call contravProjPlane2Sphere(us2(i), vs2(i), field_x_raw(i,j,iPatch2), field_y_raw(i,j,iPatch2), mesh%matrixA_ext(:,:,i,j,iPatch2))
      enddo
      us = 0.5 * ( us1 + us2(ide:ids:-1) )
      vs = 0.5 * ( vs1 + vs2(ide:ids:-1) )
      
      j = jde
      do i = ids, ide
        call contravProjSphere2Plane(field_x(i,j,iPatch1), field_y(i,j,iPatch1), us(i), vs(i), mesh%matrixIA_ext(:,:,i,j,iPatch1))
      enddo
      
      j = jde
      do i = ids, ide
        ir = ide - i + 1
        call contravProjSphere2Plane(field_x(i,j,iPatch2), field_y(i,j,iPatch2), us(ir), vs(ir), mesh%matrixIA_ext(:,:,i,j,iPatch2))
      enddo
      
      ! edge(3,6)
      iPatch1 = 3
      j = jds
      do i = ids, ide
        call contravProjPlane2Sphere(us1(i), vs1(i), field_x_raw(i,j,iPatch1), field_y_raw(i,j,iPatch1), mesh%matrixA_ext(:,:,i,j,iPatch1))
      enddo
      
      iPatch2 = 6
      j = jds
      do i = ids, ide
        call contravProjPlane2Sphere(us2(i), vs2(i), field_x_raw(i,j,iPatch2), field_y_raw(i,j,iPatch2), mesh%matrixA_ext(:,:,i,j,iPatch2))
      enddo
      us = 0.5 * ( us1 + us2(ide:ids:-1) )
      vs = 0.5 * ( vs1 + vs2(ide:ids:-1) )
      
      j = jds
      do i = ids, ide
        call contravProjSphere2Plane(field_x(i,j,iPatch1), field_y(i,j,iPatch1), us(i), vs(i), mesh%matrixIA_ext(:,:,i,j,iPatch1))
      enddo
      
      j = jds
      do i = ids, ide
        ir = ide - i + 1
        call contravProjSphere2Plane(field_x(i,j,iPatch2), field_y(i,j,iPatch2), us(ir), vs(ir), mesh%matrixIA_ext(:,:,i,j,iPatch2))
      enddo
      
      ! edge(4,5)
      iPatch1 = 4
      j = jde
      do i = ids, ide
        call contravProjPlane2Sphere(us1(i), vs1(i), field_x_raw(i,j,iPatch1), field_y_raw(i,j,iPatch1), mesh%matrixA_ext(:,:,i,j,iPatch1))
      enddo
      
      iPatch2 = 5
      i = ids
      do j = jds, jde
        call contravProjPlane2Sphere(us2(j), vs2(j), field_x_raw(i,j,iPatch2), field_y_raw(i,j,iPatch2), mesh%matrixA_ext(:,:,i,j,iPatch2))
      enddo
      us = 0.5 * ( us1 + us2(jde:jds:-1) )
      vs = 0.5 * ( vs1 + vs2(jde:jds:-1) )
      
      j = jde
      do i = ids, ide
        call contravProjSphere2Plane(field_x(i,j,iPatch1), field_y(i,j,iPatch1), us(i), vs(i), mesh%matrixIA_ext(:,:,i,j,iPatch1))
      enddo
      
      i = ids
      do j = jds, jde
        jr = jde - j + 1
        call contravProjSphere2Plane(field_x(i,j,iPatch2), field_y(i,j,iPatch2), us(jr), vs(jr), mesh%matrixIA_ext(:,:,i,j,iPatch2))
      enddo
      
      ! edge(4,6)
      iPatch1 = 4
      j = jds
      do i = ids, ide
        call contravProjPlane2Sphere(us1(i), vs1(i), field_x_raw(i,j,iPatch1), field_y_raw(i,j,iPatch1), mesh%matrixA_ext(:,:,i,j,iPatch1))
      enddo
      
      iPatch2 = 6
      i = ids
      do j = jds, jde
        call contravProjPlane2Sphere(us2(j), vs2(j), field_x_raw(i,j,iPatch2), field_y_raw(i,j,iPatch2), mesh%matrixA_ext(:,:,i,j,iPatch2))
      enddo
      us = 0.5 * ( us1 + us2 )
      vs = 0.5 * ( vs1 + vs2 )
      
      j = jds
      do i = ids, ide
        call contravProjSphere2Plane(field_x(i,j,iPatch1), field_y(i,j,iPatch1), us(i), vs(i), mesh%matrixIA_ext(:,:,i,j,iPatch1))
      enddo
      
      i = ids
      do j = jds, jde
        call contravProjSphere2Plane(field_x(i,j,iPatch2), field_y(i,j,iPatch2), us(j), vs(j), mesh%matrixIA_ext(:,:,i,j,iPatch2))
      enddo
    end subroutine unify_bdy_field_contravariant
    
    subroutine unify_bdy_field_covariant(field_x,field_y)
      real,intent(inout) :: field_x(ids:ide,jds:jde,ifs:ife)
      real,intent(inout) :: field_y(ids:ide,jds:jde,ifs:ife)
      
      real us (ids:ide),vs (ids:ide)
      real us1(ids:ide),vs1(ids:ide)
      real us2(ids:ide),vs2(ids:ide)
      
      integer i,j,iPatch
      integer ir,jr
      integer iPatch1, iPatch2
      
      real field_x_raw(ids:ide,jds:jde,ifs:ife)
      real field_y_raw(ids:ide,jds:jde,ifs:ife)
      
      field_x_raw = field_x
      field_y_raw = field_y
      
      ! Bdy lines
      ! edge(1,4)
      iPatch1 = 1
      i = ids
      do j = jds, jde
        call covProjPlane2Sphere(us1(j), vs1(j), field_x_raw(i,j,iPatch1), field_y_raw(i,j,iPatch1), mesh%matrixA_ext(:,:,i,j,iPatch1), mesh%matrixIG_ext(:,:,i,j,iPatch1))
      enddo
      
      iPatch2 = 4
      i = ide
      do j = jds, jde
        call covProjPlane2Sphere(us2(j), vs2(j), field_x_raw(i,j,iPatch2), field_y_raw(i,j,iPatch2), mesh%matrixA_ext(:,:,i,j,iPatch2), mesh%matrixIG_ext(:,:,i,j,iPatch2))
      enddo
      us = 0.5 * ( us1 + us2 )
      vs = 0.5 * ( vs1 + vs2 )
      
      i = ids
      do j = jds, jde
        call covProjSphere2Plane(field_x(i,j,iPatch1), field_y(i,j,iPatch1), us(j), vs(j), mesh%matrixIA_ext(:,:,i,j,iPatch1), mesh%matrixG_ext(:,:,i,j,iPatch1))
      enddo
      
      i = ide
      do j = jds, jde
        call covProjSphere2Plane(field_x(i,j,iPatch2), field_y(i,j,iPatch2), us(j), vs(j), mesh%matrixIA_ext(:,:,i,j,iPatch2), mesh%matrixG_ext(:,:,i,j,iPatch2))
      enddo
      

      ! edge(1,2)
      iPatch1 = 2
      i = ids
      do j = jds, jde
        call covProjPlane2Sphere(us1(j), vs1(j), field_x_raw(i,j,iPatch1), field_y_raw(i,j,iPatch1), mesh%matrixA_ext(:,:,i,j,iPatch1), mesh%matrixIG_ext(:,:,i,j,iPatch1))
      enddo
      
      iPatch2 = 1
      i = ide
      do j = jds, jde
        call covProjPlane2Sphere(us2(j), vs2(j), field_x_raw(i,j,iPatch2), field_y_raw(i,j,iPatch2), mesh%matrixA_ext(:,:,i,j,iPatch2), mesh%matrixIG_ext(:,:,i,j,iPatch2))
      enddo
      us = 0.5 * ( us1 + us2 )
      vs = 0.5 * ( vs1 + vs2 )
      
      i = ids
      do j = jds, jde
        call covProjSphere2Plane(field_x(i,j,iPatch1), field_y(i,j,iPatch1), us(j), vs(j), mesh%matrixIA_ext(:,:,i,j,iPatch1), mesh%matrixG_ext(:,:,i,j,iPatch1))
      enddo
      
      i = ide
      do j = jds, jde
        call covProjSphere2Plane(field_x(i,j,iPatch2), field_y(i,j,iPatch2), us(j), vs(j), mesh%matrixIA_ext(:,:,i,j,iPatch2), mesh%matrixG_ext(:,:,i,j,iPatch2))
      enddo
      
      ! edge(2,3)
      iPatch1 = 1
      i = ids
      do j = jds, jde
        call covProjPlane2Sphere(us1(j), vs1(j), field_x_raw(i,j,iPatch1), field_y_raw(i,j,iPatch1), mesh%matrixA_ext(:,:,i,j,iPatch1), mesh%matrixIG_ext(:,:,i,j,iPatch1))
      enddo
      
      iPatch2 = 4
      i = ide
      do j = jds, jde
        call covProjPlane2Sphere(us2(j), vs2(j), field_x_raw(i,j,iPatch2), field_y_raw(i,j,iPatch2), mesh%matrixA_ext(:,:,i,j,iPatch2), mesh%matrixIG_ext(:,:,i,j,iPatch2))
      enddo
      us = 0.5 * ( us1 + us2 )
      vs = 0.5 * ( vs1 + vs2 )
      
      i = ids
      do j = jds, jde
        call covProjSphere2Plane(field_x(i,j,iPatch1), field_y(i,j,iPatch1), us(j), vs(j), mesh%matrixIA_ext(:,:,i,j,iPatch1), mesh%matrixG_ext(:,:,i,j,iPatch1))
      enddo
      
      i = ide
      do j = jds, jde
        call covProjSphere2Plane(field_x(i,j,iPatch2), field_y(i,j,iPatch2), us(j), vs(j), mesh%matrixIA_ext(:,:,i,j,iPatch2), mesh%matrixG_ext(:,:,i,j,iPatch2))
      enddo
      
      ! edge(3,4)
      iPatch1 = 4
      i = ids
      do j = jds, jde
        call covProjPlane2Sphere(us1(j), vs1(j), field_x_raw(i,j,iPatch1), field_y_raw(i,j,iPatch1), mesh%matrixA_ext(:,:,i,j,iPatch1), mesh%matrixIG_ext(:,:,i,j,iPatch1))
      enddo
      
      iPatch2 = 3
      i = ide
      do j = jds, jde
        call covProjPlane2Sphere(us2(j), vs2(j), field_x_raw(i,j,iPatch2), field_y_raw(i,j,iPatch2), mesh%matrixA_ext(:,:,i,j,iPatch2), mesh%matrixIG_ext(:,:,i,j,iPatch2))
      enddo
      us = 0.5 * ( us1 + us2 )
      vs = 0.5 * ( vs1 + vs2 )
      
      i = ids
      do j = jds, jde
        call covProjSphere2Plane(field_x(i,j,iPatch1), field_y(i,j,iPatch1), us(j), vs(j), mesh%matrixIA_ext(:,:,i,j,iPatch1), mesh%matrixG_ext(:,:,i,j,iPatch1))
      enddo
      
      i = ide
      do j = jds, jde
        call covProjSphere2Plane(field_x(i,j,iPatch2), field_y(i,j,iPatch2), us(j), vs(j), mesh%matrixIA_ext(:,:,i,j,iPatch2), mesh%matrixG_ext(:,:,i,j,iPatch2))
      enddo
      
      ! edge(1,5)
      iPatch1 = 5
      j = jds
      do i = ids, ide
        call covProjPlane2Sphere(us1(i), vs1(i), field_x_raw(i,j,iPatch1), field_y_raw(i,j,iPatch1), mesh%matrixA_ext(:,:,i,j,iPatch1), mesh%matrixIG_ext(:,:,i,j,iPatch1))
      enddo
      
      iPatch2 = 1
      j = jde
      do i = ids, ide
        call covProjPlane2Sphere(us2(i), vs2(i), field_x_raw(i,j,iPatch2), field_y_raw(i,j,iPatch2), mesh%matrixA_ext(:,:,i,j,iPatch2), mesh%matrixIG_ext(:,:,i,j,iPatch2))
      enddo
      us = 0.5 * ( us1 + us2 )
      vs = 0.5 * ( vs1 + vs2 )
      
      j = jds
      do i = ids, ide
        call covProjSphere2Plane(field_x(i,j,iPatch1), field_y(i,j,iPatch1), us(i), vs(i), mesh%matrixIA_ext(:,:,i,j,iPatch1), mesh%matrixG_ext(:,:,i,j,iPatch1))
      enddo
      
      j = jde
      do i = ids, ide
        call covProjSphere2Plane(field_x(i,j,iPatch2), field_y(i,j,iPatch2), us(i), vs(i), mesh%matrixIA_ext(:,:,i,j,iPatch2), mesh%matrixG_ext(:,:,i,j,iPatch2))
      enddo
      
      ! edge(1,6)
      iPatch1 = 1
      j = jds
      do i = ids, ide
        call covProjPlane2Sphere(us1(i), vs1(i), field_x_raw(i,j,iPatch1), field_y_raw(i,j,iPatch1), mesh%matrixA_ext(:,:,i,j,iPatch1), mesh%matrixIG_ext(:,:,i,j,iPatch1))
      enddo
      
      iPatch2 = 6
      j = jde
      do i = ids, ide
        call covProjPlane2Sphere(us2(i), vs2(i), field_x_raw(i,j,iPatch2), field_y_raw(i,j,iPatch2), mesh%matrixA_ext(:,:,i,j,iPatch2), mesh%matrixIG_ext(:,:,i,j,iPatch2))
      enddo
      us = 0.5 * ( us1 + us2 )
      vs = 0.5 * ( vs1 + vs2 )
      
      j = jds
      do i = ids, ide
        call covProjSphere2Plane(field_x(i,j,iPatch1), field_y(i,j,iPatch1), us(i), vs(i), mesh%matrixIA_ext(:,:,i,j,iPatch1), mesh%matrixG_ext(:,:,i,j,iPatch1))
      enddo
      
      j = jde
      do i = ids, ide
        call covProjSphere2Plane(field_x(i,j,iPatch2), field_y(i,j,iPatch2), us(i), vs(i), mesh%matrixIA_ext(:,:,i,j,iPatch2), mesh%matrixG_ext(:,:,i,j,iPatch2))
      enddo
      
      ! edge(2,5)
      iPatch1 = 2
      j = jde
      do i = ids, ide
        call covProjPlane2Sphere(us1(i), vs1(i), field_x_raw(i,j,iPatch1), field_y_raw(i,j,iPatch1), mesh%matrixA_ext(:,:,i,j,iPatch1), mesh%matrixIG_ext(:,:,i,j,iPatch1))
      enddo
      
      iPatch2 = 5
      i = ide
      do j = jds, jde
        call covProjPlane2Sphere(us2(j), vs2(j), field_x_raw(i,j,iPatch2), field_y_raw(i,j,iPatch2), mesh%matrixA_ext(:,:,i,j,iPatch2), mesh%matrixIG_ext(:,:,i,j,iPatch2))
      enddo
      us = 0.5 * ( us1 + us2 )
      vs = 0.5 * ( vs1 + vs2 )
      
      j = jde
      do i = ids, ide
        call covProjSphere2Plane(field_x(i,j,iPatch1), field_y(i,j,iPatch1), us(i), vs(i), mesh%matrixIA_ext(:,:,i,j,iPatch1), mesh%matrixG_ext(:,:,i,j,iPatch1))
      enddo
      
      i = ide
      do j = jds, jde
        call covProjSphere2Plane(field_x(i,j,iPatch2), field_y(i,j,iPatch2), us(j), vs(j), mesh%matrixIA_ext(:,:,i,j,iPatch2), mesh%matrixG_ext(:,:,i,j,iPatch2))
      enddo
      
      ! edge(2,6)
      iPatch1 = 2
      j = jds
      do i = ids, ide
        call covProjPlane2Sphere(us1(i), vs1(i), field_x_raw(i,j,iPatch1), field_y_raw(i,j,iPatch1), mesh%matrixA_ext(:,:,i,j,iPatch1), mesh%matrixIG_ext(:,:,i,j,iPatch1))
      enddo
      
      iPatch2 = 6
      i = ide
      do j = jds, jde
        call covProjPlane2Sphere(us2(j), vs2(j), field_x_raw(i,j,iPatch2), field_y_raw(i,j,iPatch2), mesh%matrixA_ext(:,:,i,j,iPatch2), mesh%matrixIG_ext(:,:,i,j,iPatch2))
      enddo
      us = 0.5 * ( us1 + us2(jde:jds:-1) )
      vs = 0.5 * ( vs1 + vs2(jde:jds:-1) )
      
      j = jds
      do i = ids, ide
        call covProjSphere2Plane(field_x(i,j,iPatch1), field_y(i,j,iPatch1), us(i), vs(i), mesh%matrixIA_ext(:,:,i,j,iPatch1), mesh%matrixG_ext(:,:,i,j,iPatch1))
      enddo
      
      i = ide
      do j = jds, jde
        jr = jde - j + 1
        call covProjSphere2Plane(field_x(i,j,iPatch2), field_y(i,j,iPatch2), us(jr), vs(jr), mesh%matrixIA_ext(:,:,i,j,iPatch2), mesh%matrixG_ext(:,:,i,j,iPatch2))
      enddo
      
      ! edge(3,5)
      iPatch1 = 3
      j = jde
      do i = ids, ide
        call covProjPlane2Sphere(us1(i), vs1(i), field_x_raw(i,j,iPatch1), field_y_raw(i,j,iPatch1), mesh%matrixA_ext(:,:,i,j,iPatch1), mesh%matrixIG_ext(:,:,i,j,iPatch1))
      enddo
      
      iPatch2 = 5
      j = jde
      do i = ids, ide
        call covProjPlane2Sphere(us2(i), vs2(i), field_x_raw(i,j,iPatch2), field_y_raw(i,j,iPatch2), mesh%matrixA_ext(:,:,i,j,iPatch2), mesh%matrixIG_ext(:,:,i,j,iPatch2))
      enddo
      us = 0.5 * ( us1 + us2(ide:ids:-1) )
      vs = 0.5 * ( vs1 + vs2(ide:ids:-1) )
      
      j = jde
      do i = ids, ide
        call covProjSphere2Plane(field_x(i,j,iPatch1), field_y(i,j,iPatch1), us(i), vs(i), mesh%matrixIA_ext(:,:,i,j,iPatch1), mesh%matrixG_ext(:,:,i,j,iPatch1))
      enddo
      
      j = jde
      do i = ids, ide
        ir = ide - i + 1
        call covProjSphere2Plane(field_x(i,j,iPatch2), field_y(i,j,iPatch2), us(ir), vs(ir), mesh%matrixIA_ext(:,:,i,j,iPatch2), mesh%matrixG_ext(:,:,i,j,iPatch2))
      enddo
      
      ! edge(3,6)
      iPatch1 = 3
      j = jds
      do i = ids, ide
        call covProjPlane2Sphere(us1(i), vs1(i), field_x_raw(i,j,iPatch1), field_y_raw(i,j,iPatch1), mesh%matrixA_ext(:,:,i,j,iPatch1), mesh%matrixIG_ext(:,:,i,j,iPatch1))
      enddo
      
      iPatch2 = 6
      j = jds
      do i = ids, ide
        call covProjPlane2Sphere(us2(i), vs2(i), field_x_raw(i,j,iPatch2), field_y_raw(i,j,iPatch2), mesh%matrixA_ext(:,:,i,j,iPatch2), mesh%matrixIG_ext(:,:,i,j,iPatch2))
      enddo
      us = 0.5 * ( us1 + us2(ide:ids:-1) )
      vs = 0.5 * ( vs1 + vs2(ide:ids:-1) )
      
      j = jds
      do i = ids, ide
        call covProjSphere2Plane(field_x(i,j,iPatch1), field_y(i,j,iPatch1), us(i), vs(i), mesh%matrixIA_ext(:,:,i,j,iPatch1), mesh%matrixG_ext(:,:,i,j,iPatch1))
      enddo
      
      j = jds
      do i = ids, ide
        ir = ide - i + 1
        call covProjSphere2Plane(field_x(i,j,iPatch2), field_y(i,j,iPatch2), us(ir), vs(ir), mesh%matrixIA_ext(:,:,i,j,iPatch2), mesh%matrixG_ext(:,:,i,j,iPatch2))
      enddo
      
      ! edge(4,5)
      iPatch1 = 4
      j = jde
      do i = ids, ide
        call covProjPlane2Sphere(us1(i), vs1(i), field_x_raw(i,j,iPatch1), field_y_raw(i,j,iPatch1), mesh%matrixA_ext(:,:,i,j,iPatch1), mesh%matrixIG_ext(:,:,i,j,iPatch1))
      enddo
      
      iPatch2 = 5
      i = ids
      do j = jds, jde
        call covProjPlane2Sphere(us2(j), vs2(j), field_x_raw(i,j,iPatch2), field_y_raw(i,j,iPatch2), mesh%matrixA_ext(:,:,i,j,iPatch2), mesh%matrixIG_ext(:,:,i,j,iPatch2))
      enddo
      us = 0.5 * ( us1 + us2(jde:jds:-1) )
      vs = 0.5 * ( vs1 + vs2(jde:jds:-1) )
      
      j = jde
      do i = ids, ide
        call covProjSphere2Plane(field_x(i,j,iPatch1), field_y(i,j,iPatch1), us(i), vs(i), mesh%matrixIA_ext(:,:,i,j,iPatch1), mesh%matrixG_ext(:,:,i,j,iPatch1))
      enddo
      
      i = ids
      do j = jds, jde
        jr = jde - j + 1
        call covProjSphere2Plane(field_x(i,j,iPatch2), field_y(i,j,iPatch2), us(jr), vs(jr), mesh%matrixIA_ext(:,:,i,j,iPatch2), mesh%matrixG_ext(:,:,i,j,iPatch2))
      enddo
      
      ! edge(4,6)
      iPatch1 = 4
      j = jds
      do i = ids, ide
        call covProjPlane2Sphere(us1(i), vs1(i), field_x_raw(i,j,iPatch1), field_y_raw(i,j,iPatch1), mesh%matrixA_ext(:,:,i,j,iPatch1), mesh%matrixIG_ext(:,:,i,j,iPatch1))
      enddo
      
      iPatch2 = 6
      i = ids
      do j = jds, jde
        call covProjPlane2Sphere(us2(j), vs2(j), field_x_raw(i,j,iPatch2), field_y_raw(i,j,iPatch2), mesh%matrixA_ext(:,:,i,j,iPatch2), mesh%matrixIG_ext(:,:,i,j,iPatch2))
      enddo
      us = 0.5 * ( us1 + us2 )
      vs = 0.5 * ( vs1 + vs2 )
      
      j = jds
      do i = ids, ide
        call covProjSphere2Plane(field_x(i,j,iPatch1), field_y(i,j,iPatch1), us(i), vs(i), mesh%matrixIA_ext(:,:,i,j,iPatch1), mesh%matrixG_ext(:,:,i,j,iPatch1))
      enddo
      
      i = ids
      do j = jds, jde
        call covProjSphere2Plane(field_x(i,j,iPatch2), field_y(i,j,iPatch2), us(j), vs(j), mesh%matrixIA_ext(:,:,i,j,iPatch2), mesh%matrixG_ext(:,:,i,j,iPatch2))
      enddo
    end subroutine unify_bdy_field_covariant
    
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

