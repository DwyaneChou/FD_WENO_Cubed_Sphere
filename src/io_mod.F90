module io_mod
  use netcdf
  use constants_mod
  use parameters_mod
  use mesh_mod
  use stat_mod
  use tend_mod
  use math_mod
  use ghost_mod
  implicit none
    
    character(13) :: ncFile = 'mcv_output.nc'
    
    contains
    subroutine read_mesh
      
      integer ncid,status
      
      integer lon_id,lat_id
      integer x_id,y_id,z_id
      integer xi_id,eta_id
      integer jab_id
      integer sqrtG_id
      integer matrixG_id,matrixIG_id
      integer matrixA_id,matrixIA_id
      
      integer iCell,jCell,iPatch
      integer ic,jc
      integer i,j,k
      integer stat
      
      real lon,lat
      real dlondx,dlatdx
      real dlondy,dlatdy
      real dlondz,dlatdz
      real sph2cart_matrix(2,3)
      real IA1(2,2)
      real IA2(2,2)
      
      status = nf90_open(trim(mesh_file),NF90_NOWRITE,ncid)
      call handle_err(status)
      
      status = nf90_get_att(ncid,NF90_GLOBAL,'dx'        ,dx        )
      status = nf90_get_att(ncid,NF90_GLOBAL,'dy'        ,dy        )
      status = nf90_get_att(ncid,NF90_GLOBAL,'ids'       ,ids       )
      status = nf90_get_att(ncid,NF90_GLOBAL,'ide'       ,ide       )
      status = nf90_get_att(ncid,NF90_GLOBAL,'jds'       ,jds       )
      status = nf90_get_att(ncid,NF90_GLOBAL,'jde'       ,jde       )
      status = nf90_get_att(ncid,NF90_GLOBAL,'ifs'       ,ifs       )
      status = nf90_get_att(ncid,NF90_GLOBAL,'ife'       ,ife       )
      
      its = ids
      ite = ( ide - ids ) / 2
      jts = jds
      jte = ( jde - jds ) / 2
      
      dx = dx * D2R
      dy = dy * D2R
      
      ! Calculate element numbers on x/y direction
      Nx = ite
      Ny = jte
      
      ! Calculate PV number on x/y direction
      nPVx = Nx * (DOF - 1) + 1
      nPVy = Ny * (DOF - 1) + 1
      
      ! Calculate grid number on longitude/latitude coordinate
      Nlambda = nPVx
      Ntheta  = nPVy
      
      ! Calculate starting and ending index for physical domain
      ics  = 1  - xhalo
      ice  = Nx + xhalo
      jcs  = 1  - yhalo
      jce  = Ny + yhalo
      
      Nx_halo = ice - ics + 1
      Ny_halo = jce - jcs + 1
      
      ! Calculate starting and ending index for memory array
      ips  = 1    - xhalo * (DOF - 1)
      ipe  = nPVx + xhalo * (DOF - 1)
      jps  = 1    - yhalo * (DOF - 1)
      jpe  = nPVy + yhalo * (DOF - 1)
      
      nPVx_halo = ipe - ips + 1
      nPVy_halo = jpe - jps + 1
      
      nPVHalo = xhalo * (DOF - 1)
      
      call initMesh
      
      status = nf90_inq_varid(ncid,'lon'        ,lon_id     )
      status = nf90_inq_varid(ncid,'lat'        ,lat_id     )
      status = nf90_inq_varid(ncid,'xi'         ,xi_id      )
      status = nf90_inq_varid(ncid,'eta'        ,eta_id     )
      status = nf90_inq_varid(ncid,'x'          ,x_id       )
      status = nf90_inq_varid(ncid,'y'          ,y_id       )
      status = nf90_inq_varid(ncid,'z'          ,z_id       )
      status = nf90_inq_varid(ncid,'jab'        ,jab_id     )
      
      status = nf90_get_var(ncid,lon_id     , mesh%lon_ext     )
      status = nf90_get_var(ncid,lat_id     , mesh%lat_ext     )
      status = nf90_get_var(ncid,xi_id      , mesh%xi_ext      )
      status = nf90_get_var(ncid,eta_id     , mesh%eta_ext     )
      status = nf90_get_var(ncid,x_id       , mesh%x_ext       )
      status = nf90_get_var(ncid,y_id       , mesh%y_ext       )
      status = nf90_get_var(ncid,z_id       , mesh%z_ext       )
      status = nf90_get_var(ncid,jab_id     , mesh%jab_ext     )
      if(status/=nf90_noerr) call handle_err(status)
      
      mesh%lon_ext = mesh%lon_ext * D2R
      mesh%lat_ext = mesh%lat_ext * D2R
      
      mesh%x_ext = mesh%x_ext * radius
      mesh%y_ext = mesh%y_ext * radius
      mesh%z_ext = mesh%z_ext * radius
      
      mesh%jab_ext      = mesh%jab_ext   * radius
      
      mesh%lon = FillValue
      mesh%lat = FillValue
      do iPatch = ifs,ife
        do jCell = jts,jte
          do iCell = its,ite
            ic = iCell * 2
            jc = jCell * 2
            mesh%lon     (    iCell,jCell,iPatch) = mesh%lon_ext     (    ic,jc,iPatch)
            mesh%lat     (    iCell,jCell,iPatch) = mesh%lat_ext     (    ic,jc,iPatch)
            mesh%xi      (    iCell,jCell,iPatch) = mesh%xi_ext      (    ic,jc,iPatch)
            mesh%eta     (    iCell,jCell,iPatch) = mesh%eta_ext     (    ic,jc,iPatch)
            mesh%x       (    iCell,jCell,iPatch) = mesh%x_ext       (    ic,jc,iPatch)
            mesh%y       (    iCell,jCell,iPatch) = mesh%y_ext       (    ic,jc,iPatch)
            mesh%z       (    iCell,jCell,iPatch) = mesh%z_ext       (    ic,jc,iPatch)
            mesh%jab     (:,:,iCell,jCell,iPatch) = mesh%jab_ext     (:,:,ic,jc,iPatch)
          enddo
        enddo
      enddo
      
      mesh%sinlon = sin(mesh%lon)
      mesh%coslon = cos(mesh%lon)
      mesh%sinlat = sin(mesh%lat)
      mesh%coslat = cos(mesh%lat)
      
      call CubedSphereFillHalo(mesh%lon   )
      call CubedSphereFillHalo(mesh%lat   )
      call CubedSphereFillHalo(mesh%x     )
      call CubedSphereFillHalo(mesh%y     )
      call CubedSphereFillHalo(mesh%z     )
      call CubedSphereFillHalo(mesh%sinlon)
      call CubedSphereFillHalo(mesh%coslon)
      call CubedSphereFillHalo(mesh%sinlat)
      call CubedSphereFillHalo(mesh%coslat)
      
      call CubedSphereFillJab(mesh%jab)
      
      ! Reset matrices
      do k = ifs,ife
        do j = jcs,jce
          do i = ics,ice
            lon = mesh%lon(i,j,k)
            lat = mesh%lat(i,j,k)
            
            if(lon/=FillValue .and. lat/=FillValue)then
              dlondx = -sin(lon); dlatdx = -cos(lon)*sin(lat)
              dlondy =  cos(lon); dlatdy = -sin(lon)*sin(lat)
              dlondz = 0.       ; dlatdz =  cos(lat)
              
              sph2cart_matrix(1,1) = dlondx
              sph2cart_matrix(1,2) = dlondy
              sph2cart_matrix(1,3) = dlondz
              sph2cart_matrix(2,1) = dlatdx
              sph2cart_matrix(2,2) = dlatdy
              sph2cart_matrix(2,3) = dlatdz
              
              mesh%matrixG(:,:,i,j,k) = matmul( transpose(mesh%jab(:,:,i,j,k)), mesh%jab(:,:,i,j,k) )
              
              call BRINV(2,mesh%matrixG(:,:,i,j,k),mesh%matrixIG(:,:,i,j,k),stat)
              if(stat==0)stop 'Fail to calculate inverse of matrixG'
              
              mesh%sqrtG(i,j,k) = sqrt(det(mesh%matrixG(:,:,i,j,k)))
              
              mesh%matrixA(:,:,i,j,k) = matmul(sph2cart_matrix,mesh%jab(:,:,i,j,k))
              
              if(abs(lat*R2D)/=90)then
                call BRINV(2,mesh%matrixA(:,:,i,j,k),mesh%matrixIA(:,:,i,j,k),stat)
                if(stat==0)then
                  print*,'Fail to calculate inverse of matrixA at i,j,k,det :',i,j,k,det(mesh%matrixA(:,:,i,j,k))
                  stop
                endif
              elseif(lat*R2D==90)then
                ! Replace by analytical approximation on north pole
                IA1(1,1) = cos(lon)
                IA1(1,2) = -sin(lon)
                IA1(2,1) = sin(lon)
                IA1(2,2) = cos(lon)
                
                IA2(1,1) = 1.
                IA2(1,2) = 0.
                IA2(2,1) = 0.
                IA2(2,2) = 1.
                
                mesh%matrixIA(:,:,i,j,k) = matmul(IA2,IA1) * radius
                !print*,mesh%matrixIA(:,:,i,j,k)
              elseif(lat*R2D==-90)then
                ! Replace by analytical approximation on south pole
                IA1(1,1) = -cos(lon)
                IA1(1,2) = sin(lon)
                IA1(2,1) = sin(lon)
                IA1(2,2) = cos(lon)
                
                IA2(1,1) = 1.
                IA2(1,2) = 0.
                IA2(2,1) = 0.
                IA2(2,2) = 1.
                
                mesh%matrixIA(:,:,i,j,k) = matmul(IA2,IA1) * radius
              endif
            endif
          enddo
        enddo
      enddo
      
      do k = ifs,ife
        do j = jds,jde
          do i = ids,ide
            lon = mesh%lon_ext(i,j,k)
            lat = mesh%lat_ext(i,j,k)
            
            if(lon/=FillValue .and. lat/=FillValue)then
              dlondx = -sin(lon); dlatdx = -cos(lon)*sin(lat)
              dlondy =  cos(lon); dlatdy = -sin(lon)*sin(lat)
              dlondz = 0.       ; dlatdz =  cos(lat)
              
              sph2cart_matrix(1,1) = dlondx
              sph2cart_matrix(1,2) = dlondy
              sph2cart_matrix(1,3) = dlondz
              sph2cart_matrix(2,1) = dlatdx
              sph2cart_matrix(2,2) = dlatdy
              sph2cart_matrix(2,3) = dlatdz
              
              mesh%matrixG_ext(:,:,i,j,k) = matmul( transpose(mesh%jab_ext(:,:,i,j,k)), mesh%jab_ext(:,:,i,j,k) )
              
              call BRINV(2,mesh%matrixG_ext(:,:,i,j,k),mesh%matrixIG_ext(:,:,i,j,k),stat)
              if(stat==0)stop 'Fail to calculate inverse of matrixG'
              
              mesh%sqrtG_ext(i,j,k) = sqrt(det(mesh%matrixG_ext(:,:,i,j,k)))
              
              mesh%matrixA_ext(:,:,i,j,k) = matmul(sph2cart_matrix,mesh%jab_ext(:,:,i,j,k))
              
              if(abs(lat*R2D)/=90)then
                call BRINV(2,mesh%matrixA_ext(:,:,i,j,k),mesh%matrixIA_ext(:,:,i,j,k),stat)
                if(stat==0)then
                  print*,'Fail to calculate inverse of matrixA at i,j,k,det :',i,j,k,det(mesh%matrixA_ext(:,:,i,j,k))
                  stop
                endif
              elseif(lat*R2D==90)then
                ! Replace by analytical approximation on north pole
                IA1(1,1) = cos(lon)
                IA1(1,2) = -sin(lon)
                IA1(2,1) = sin(lon)
                IA1(2,2) = cos(lon)
                
                IA2(1,1) = 1.
                IA2(1,2) = 0.
                IA2(2,1) = 0.
                IA2(2,2) = 1.
                
                mesh%matrixIA_ext(:,:,i,j,k) = matmul(IA2,IA1) * radius
                !print*,mesh%matrixIA(:,:,i,j,k)
              elseif(lat*R2D==-90)then
                ! Replace by analytical approximation on south pole
                IA1(1,1) = -cos(lon)
                IA1(1,2) = sin(lon)
                IA1(2,1) = sin(lon)
                IA1(2,2) = cos(lon)
                
                IA2(1,1) = 1.
                IA2(1,2) = 0.
                IA2(2,1) = 0.
                IA2(2,2) = 1.
                
                mesh%matrixIA_ext(:,:,i,j,k) = matmul(IA2,IA1) * radius
              endif
            endif
          enddo
        enddo
      enddo
      
      ! Calculate mesh infomation on VIA
      do iPatch = ifs, ife
        do jCell = jts, jte
          do iCell = its, ite
            mesh%f(iCell, jCell, iPatch) = 2. * Omega * mesh%sinlat(iCell, jCell, iPatch)
          end do
        end do
      end do
      
      ! Calculate areaCell
      mesh%areaCell = mesh%sqrtG * dx**2
      
      print*,'Spherical area, numerical, analytical :',sum(mesh%areaCell(its:ite,jts:jte,:)),4.*pi*radius**2
    end subroutine read_mesh
    
    subroutine history_init(stat)
      type(stat_field), intent(in) :: stat
      
      integer status
      integer ncid
      integer lon_dim_id,lat_dim_id
      integer patch_dim_id
      integer time_dim_id
      integer lon_id,lat_id
      integer x_id,y_id
      integer areaCell_id
      integer time_id
      integer u_id,v_id
      integer uc_id,vc_id
      integer phi_id
      integer phis_id
      integer phit_id
      integer zonal_wind_id,meridional_wind_id
      
      status = nf90_create(ncFile, NF90_CLOBBER + NF90_64BIT_OFFSET , ncid)
      if(status/=nf90_noerr) call handle_err(status)
      
      status = nf90_def_dim(ncid,'lon'   ,Nx_halo       ,lon_dim_id )
      status = nf90_def_dim(ncid,'lat'   ,Ny_halo       ,lat_dim_id )
      status = nf90_def_dim(ncid,'nPatch',Nf            ,patch_dim_id)
      status = nf90_def_dim(ncid,'time'  ,NF90_UNLIMITED,time_dim_id )
      if(status/=nf90_noerr) call handle_err(status)
      
      status = nf90_def_var(ncid,'lon'            ,NF90_DOUBLE,(/lon_dim_id,lat_dim_id,patch_dim_id            /),lon_id            )
      status = nf90_def_var(ncid,'lat'            ,NF90_DOUBLE,(/lon_dim_id,lat_dim_id,patch_dim_id            /),lat_id            )
      status = nf90_def_var(ncid,'areaCell'       ,NF90_DOUBLE,(/lon_dim_id,lat_dim_id,patch_dim_id            /),areaCell_id       )
      status = nf90_def_var(ncid,'phis'           ,NF90_DOUBLE,(/lon_dim_id,lat_dim_id,patch_dim_id            /),phis_id           )
      status = nf90_def_var(ncid,'time'           ,NF90_INT   ,(/                                   time_dim_id/),time_id           )
      status = nf90_def_var(ncid,'phi'            ,NF90_DOUBLE,(/lon_dim_id,lat_dim_id,patch_dim_id,time_dim_id/),phi_id            )
      status = nf90_def_var(ncid,'phit'           ,NF90_DOUBLE,(/lon_dim_id,lat_dim_id,patch_dim_id,time_dim_id/),phit_id           )
      status = nf90_def_var(ncid,'u'              ,NF90_DOUBLE,(/lon_dim_id,lat_dim_id,patch_dim_id,time_dim_id/),u_id              )
      status = nf90_def_var(ncid,'v'              ,NF90_DOUBLE,(/lon_dim_id,lat_dim_id,patch_dim_id,time_dim_id/),v_id              )
      status = nf90_def_var(ncid,'uc'             ,NF90_DOUBLE,(/lon_dim_id,lat_dim_id,patch_dim_id,time_dim_id/),uc_id             )
      status = nf90_def_var(ncid,'vc'             ,NF90_DOUBLE,(/lon_dim_id,lat_dim_id,patch_dim_id,time_dim_id/),vc_id             )
      status = nf90_def_var(ncid,'zonal_wind'     ,NF90_DOUBLE,(/lon_dim_id,lat_dim_id,patch_dim_id,time_dim_id/),zonal_wind_id     )
      status = nf90_def_var(ncid,'meridional_wind',NF90_DOUBLE,(/lon_dim_id,lat_dim_id,patch_dim_id,time_dim_id/),meridional_wind_id)
      if(status/=nf90_noerr) call handle_err(status)
      
      !print*,'nf90_put_att'
      status = nf90_put_att(ncid,nf90_global       ,'MCV_ORDER',DOF)
      status = nf90_put_att(ncid,nf90_global       ,'dx'       ,dx*R2D)
      status = nf90_put_att(ncid,nf90_global       ,'dy'       ,dy*R2D)
      status = nf90_put_att(ncid,nf90_global       ,'dt'       ,dt)
      status = nf90_put_att(ncid,nf90_global       ,'xhalo'    ,xhalo)
      status = nf90_put_att(ncid,nf90_global       ,'yhalo'    ,yhalo)
      status = nf90_put_att(ncid,nf90_global       ,'case_num' ,case_num)
      status = nf90_put_att(ncid,nf90_global       ,'ics'      ,ics)
      status = nf90_put_att(ncid,nf90_global       ,'ice'      ,ice)
      status = nf90_put_att(ncid,nf90_global       ,'jcs'      ,jcs)
      status = nf90_put_att(ncid,nf90_global       ,'jce'      ,jce)
      status = nf90_put_att(ncid,nf90_global       ,'ips'      ,ips)
      status = nf90_put_att(ncid,nf90_global       ,'ipe'      ,ipe)
      status = nf90_put_att(ncid,nf90_global       ,'jps'      ,jps)
      status = nf90_put_att(ncid,nf90_global       ,'jpe'      ,jpe)
      status = nf90_put_att(ncid,nf90_global       ,'its'      ,its)
      status = nf90_put_att(ncid,nf90_global       ,'ite'      ,ite)
      status = nf90_put_att(ncid,nf90_global       ,'jts'      ,jts)
      status = nf90_put_att(ncid,nf90_global       ,'jte'      ,jte)
      status = nf90_put_att(ncid,nf90_global       ,'ifs'      ,ifs)
      status = nf90_put_att(ncid,nf90_global       ,'ife'      ,ife)
      
      status = nf90_put_att(ncid,lon_id            ,'units'    ,'degree_east' )
      status = nf90_put_att(ncid,lat_id            ,'units'    ,'degree_north')
      status = nf90_put_att(ncid,areaCell_id       ,'units'    ,'m^2'         )
      status = nf90_put_att(ncid,time_id           ,'units'    ,'seconds'     )
      status = nf90_put_att(ncid,u_id              ,'units'    ,'m/s'         )
      status = nf90_put_att(ncid,v_id              ,'units'    ,'m/s'         )
      status = nf90_put_att(ncid,uc_id             ,'units'    ,'m/s'         )
      status = nf90_put_att(ncid,vc_id             ,'units'    ,'m/s'         )
      status = nf90_put_att(ncid,phi_id            ,'units'    ,'m^2/s^2'     )
      status = nf90_put_att(ncid,phis_id           ,'units'    ,'m^2/s^2'     )
      status = nf90_put_att(ncid,phit_id           ,'units'    ,'m^2/s^2'     )
      status = nf90_put_att(ncid,zonal_wind_id     ,'units'    ,'m/s'         )
      status = nf90_put_att(ncid,meridional_wind_id,'units'    ,'m/s'         )
      
      status = nf90_put_att(ncid,lon_id            ,'long_name','longitude on sphere coordinate for Cells' )
      status = nf90_put_att(ncid,lat_id            ,'long_name','latitude on sphere coordinate for Cells'  )
      status = nf90_put_att(ncid,areaCell_id       ,'long_name','area of cells'                            )
      status = nf90_put_att(ncid,time_id           ,'long_name','time'                                     )
      status = nf90_put_att(ncid,u_id              ,'long_name','covariant wind vector on x direcvtion'    )
      status = nf90_put_att(ncid,v_id              ,'long_name','covariant wind vector on y direcvtion'    )
      status = nf90_put_att(ncid,uc_id             ,'long_name','contravariant wind vector on x direcvtion')
      status = nf90_put_att(ncid,vc_id             ,'long_name','contravariant wind vector on y direcvtion')
      status = nf90_put_att(ncid,phi_id            ,'long_name','geopotential height on points'            )
      status = nf90_put_att(ncid,phis_id           ,'long_name','surface height on points'                 )
      status = nf90_put_att(ncid,phit_id           ,'long_name','total geopotential height on points'      )
      status = nf90_put_att(ncid,zonal_wind_id     ,'long_name','zonal wind'                               )
      status = nf90_put_att(ncid,meridional_wind_id,'long_name','meridional wind'                          )
      
      ! Define coordinates
      status = nf90_put_att(ncid, x_id              ,'_CoordinateAxisTypes','lon lat nPatch')
      status = nf90_put_att(ncid, y_id              ,'_CoordinateAxisTypes','lon lat nPatch')
      status = nf90_put_att(ncid, lon_id            ,'_CoordinateAxisTypes','lon lat nPatch')
      status = nf90_put_att(ncid, lat_id            ,'_CoordinateAxisTypes','lon lat nPatch')
      status = nf90_put_att(ncid, phis_id           ,'_CoordinateAxisTypes','lon lat nPatch')
      status = nf90_put_att(ncid, u_id              ,'_CoordinateAxisTypes','lon lat nPatch time')
      status = nf90_put_att(ncid, v_id              ,'_CoordinateAxisTypes','lon lat nPatch time')
      status = nf90_put_att(ncid, uc_id             ,'_CoordinateAxisTypes','lon lat nPatch time')
      status = nf90_put_att(ncid, vc_id             ,'_CoordinateAxisTypes','lon lat nPatch time')
      status = nf90_put_att(ncid, phi_id            ,'_CoordinateAxisTypes','lon lat nPatch time')
      status = nf90_put_att(ncid, phit_id           ,'_CoordinateAxisTypes','lon lat nPatch time')
      status = nf90_put_att(ncid, areaCell_id       ,'_CoordinateAxisTypes','lon lat nPatch time')
      status = nf90_put_att(ncid, zonal_wind_id     ,'_CoordinateAxisTypes','lon lat nPatch time')
      status = nf90_put_att(ncid, meridional_wind_id,'_CoordinateAxisTypes','lon lat nPatch time')
      
      status = nf90_put_att(ncid,u_id              ,'_FillValue',FillValue)
      status = nf90_put_att(ncid,v_id              ,'_FillValue',FillValue)
      status = nf90_put_att(ncid,uc_id             ,'_FillValue',FillValue)
      status = nf90_put_att(ncid,vc_id             ,'_FillValue',FillValue)
      status = nf90_put_att(ncid,phi_id            ,'_FillValue',FillValue)
      status = nf90_put_att(ncid,phis_id           ,'_FillValue',FillValue)
      status = nf90_put_att(ncid,phit_id           ,'_FillValue',FillValue)
      status = nf90_put_att(ncid,zonal_wind_id     ,'_FillValue',FillValue)
      status = nf90_put_att(ncid,meridional_wind_id,'_FillValue',FillValue)
      if(status/=nf90_noerr) call handle_err(status)
      
      status = nf90_enddef(ncid)
      if(status/=nf90_noerr) call handle_err(status)
      
      status = nf90_put_var(ncid,lon_id     , mesh%lon      * R2D)
      status = nf90_put_var(ncid,lat_id     , mesh%lat      * R2D)
      status = nf90_put_var(ncid,areaCell_id, mesh%areaCell      )
      status = nf90_put_var(ncid,phis_id    , mesh%phis          )
      !status = nf90_put_var(ncid,lon_id     , mesh%lon     (its:ite,jts:jte,ifs:ife) * R2D)
      !status = nf90_put_var(ncid,lat_id     , mesh%lat     (its:ite,jts:jte,ifs:ife) * R2D)
      !status = nf90_put_var(ncid,areaCell_id, mesh%areaCell(its:ite,jts:jte,ifs:ife)      )
      !status = nf90_put_var(ncid,phis_id    , mesh%phis    (its:ite,jts:jte,ifs:ife)      )
      if(status/=nf90_noerr) call handle_err(status)
      
      status = nf90_close(ncid)
      if(status/=nf90_noerr) call handle_err(status)
      
    end subroutine history_init
    
    subroutine history_write_stat(stat,time_slot_num)
      type(stat_field), intent(in) :: stat
      integer         , intent(in) :: time_slot_num
      
      integer status
      integer ncid
      integer time_id
      integer u_id,v_id
      integer uc_id,vc_id
      integer phi_id
      integer phit_id
      integer zonal_wind_id,meridional_wind_id
      
      !real, dimension(its:ite,jts:jte,ifs:ife) :: varout
      real, dimension(ics:ice,jcs:jce,ifs:ife) :: varout
      
      integer :: time(1)
      
      integer i,j,iPatch
      
      time(1) = time_slot_num
      !print*,'nf90_open'
      status = nf90_open(ncFile,NF90_WRITE,ncid)
      if(status/=nf90_noerr) call handle_err(status)
      
      !print*,'nf90_inq_varid'
      status = nf90_inq_varid(ncid,'time'           , time_id           )
      status = nf90_inq_varid(ncid,'u'              , u_id              )
      status = nf90_inq_varid(ncid,'v'              , v_id              )
      status = nf90_inq_varid(ncid,'uc'             , uc_id             )
      status = nf90_inq_varid(ncid,'vc'             , vc_id             )
      status = nf90_inq_varid(ncid,'phi'            , phi_id            )
      status = nf90_inq_varid(ncid,'phit'           , phit_id           )
      status = nf90_inq_varid(ncid,'zonal_wind'     , zonal_wind_id     )
      status = nf90_inq_varid(ncid,'meridional_wind', meridional_wind_id)
      if(status/=nf90_noerr) call handle_err(status)
      
      !print*,'nf90_put_var'
      status = nf90_put_var(ncid, time_id, time, start=(/time_slot_num/),count=(/1/))
      
      ! phi
      varout = stat%phi
      status = nf90_put_var(ncid, phi_id, varout, start=(/1,1,1,time_slot_num/),count=(/Nx_halo,Ny_halo,Nf,1/))
      !varout = stat%phi(its:ite,jts:jte,ifs:ife)
      !status = nf90_put_var(ncid, phi_id, varout, start=(/1,1,1,time_slot_num/),count=(/Nx,Ny,Nf,1/))
      
      !u
      varout = stat%u
      status = nf90_put_var(ncid, u_id, varout, start=(/1,1,1,time_slot_num/),count=(/Nx_halo,Ny_halo,Nf,1/))
      !varout = stat%u(its:ite,jts:jte,ifs:ife)
      !status = nf90_put_var(ncid, u_id, varout, start=(/1,1,1,time_slot_num/),count=(/Nx,Ny,Nf,1/))
      
      !v
      varout = stat%v
      status = nf90_put_var(ncid, v_id, varout, start=(/1,1,1,time_slot_num/),count=(/Nx_halo,Ny_halo,Nf,1/))
      !varout = stat%v(its:ite,jts:jte,ifs:ife)
      !status = nf90_put_var(ncid, v_id, varout, start=(/1,1,1,time_slot_num/),count=(/Nx,Ny,Nf,1/))
      
      !contra u
      varout = stat%uC
      status = nf90_put_var(ncid, uc_id, varout, start=(/1,1,1,time_slot_num/),count=(/Nx_halo,Ny_halo,Nf,1/))
      !varout = stat%uC(its:ite,jts:jte,ifs:ife)
      !status = nf90_put_var(ncid, uc_id, varout, start=(/1,1,1,time_slot_num/),count=(/Nx,Ny,Nf,1/))
      
      !contra v
      varout = stat%vC
      status = nf90_put_var(ncid, vc_id, varout, start=(/1,1,1,time_slot_num/),count=(/Nx_halo,Ny_halo,Nf,1/))
      !varout = stat%vC(its:ite,jts:jte,ifs:ife)
      !status = nf90_put_var(ncid, vc_id, varout, start=(/1,1,1,time_slot_num/),count=(/Nx,Ny,Nf,1/))
      
      !phit
      varout = stat%phi + mesh%phis
      status = nf90_put_var(ncid,phit_id, varout , start=(/1,1,1,time_slot_num/),count=(/Nx_halo,Ny_halo,Nf,1/))
      !varout = stat%phi(its:ite,jts:jte,ifs:ife) + mesh%phis(its:ite,jts:jte,ifs:ife)
      !status = nf90_put_var(ncid,phit_id, varout , start=(/1,1,1,time_slot_num/),count=(/Nx,Ny,Nf,1/))
      
      ! zonal_wind
      varout = stat%zonal_wind
      status = nf90_put_var(ncid,zonal_wind_id, varout, start=(/1,1,1,time_slot_num/), count=(/Nx_halo,Ny_halo,Nf,1/))
      !varout = stat%zonal_wind(its:ite,jts:jte,ifs:ife)
      !status = nf90_put_var(ncid,zonal_wind_id, varout, start=(/1,1,1,time_slot_num/), count=(/Nx,Ny,Nf,1/))
      
      ! meridional_wind
      varout = stat%meridional_wind
      status = nf90_put_var(ncid,meridional_wind_id, varout, start=(/1,1,1,time_slot_num/), count=(/Nx_halo,Ny_halo,Nf,1/))
      !varout = stat%meridional_wind(its:ite,jts:jte,ifs:ife)
      !status = nf90_put_var(ncid,meridional_wind_id, varout, start=(/1,1,1,time_slot_num/), count=(/Nx,Ny,Nf,1/))
      if(status/=nf90_noerr) call handle_err(status)
      
      !print*,'nf90_close'
      status = nf90_close(ncid)
      if(status/=nf90_noerr) call handle_err(status)
      
    end subroutine history_write_stat
    
    subroutine history_write_tend(tend,time_slot_num)
      type(tend_field), intent(in) :: tend
      integer         , intent(in) :: time_slot_num
      
      integer status
      integer ncid
      integer time_id
      integer u_id,v_id
      integer uc_id,vc_id
      integer phi_id
      integer phit_id
      integer zonal_wind_id,meridional_wind_id
      
      !real, dimension(its:ite,jts:jte,ifs:ife) :: varout
      real, dimension(ics:ice,jcs:jce,ifs:ife) :: varout
      
      integer :: time(1)
      
      integer i,j,iPatch
      
      time(1) = time_slot_num
      !print*,'nf90_open'
      status = nf90_open(ncFile,NF90_WRITE,ncid)
      if(status/=nf90_noerr) call handle_err(status)
      
      !print*,'nf90_inq_varid'
      status = nf90_inq_varid(ncid,'time'           , time_id           )
      status = nf90_inq_varid(ncid,'u'              , u_id              )
      status = nf90_inq_varid(ncid,'v'              , v_id              )
      status = nf90_inq_varid(ncid,'phit'           , phit_id           )
      if(status/=nf90_noerr) call handle_err(status)
      
      !print*,'nf90_put_var'
      status = nf90_put_var(ncid, time_id, time, start=(/time_slot_num/),count=(/1/))
      
      ! phi
      varout = tend%phiG
      status = nf90_put_var(ncid, phi_id, varout, start=(/1,1,1,time_slot_num/),count=(/Nx_halo,Ny_halo,Nf,1/))
      !varout = stat%phi(its:ite,jts:jte,ifs:ife)
      !status = nf90_put_var(ncid, phi_id, varout, start=(/1,1,1,time_slot_num/),count=(/Nx,Ny,Nf,1/))
      
      !u
      varout = tend%u
      status = nf90_put_var(ncid, u_id, varout, start=(/1,1,1,time_slot_num/),count=(/Nx_halo,Ny_halo,Nf,1/))
      !varout = stat%u(its:ite,jts:jte,ifs:ife)
      !status = nf90_put_var(ncid, u_id, varout, start=(/1,1,1,time_slot_num/),count=(/Nx,Ny,Nf,1/))
      
      !v
      varout = tend%v
      status = nf90_put_var(ncid, v_id, varout, start=(/1,1,1,time_slot_num/),count=(/Nx_halo,Ny_halo,Nf,1/))
      !varout = stat%v(its:ite,jts:jte,ifs:ife)
      !status = nf90_put_var(ncid, v_id, varout, start=(/1,1,1,time_slot_num/),count=(/Nx,Ny,Nf,1/))
      
      if(status/=nf90_noerr) call handle_err(status)
      
      !print*,'nf90_close'
      status = nf90_close(ncid)
      if(status/=nf90_noerr) call handle_err(status)
      
    end subroutine history_write_tend
    
    subroutine handle_err(status)
      implicit none
      integer,intent(in)::status
            
      if(status/=nf90_noerr)then
          print*, trim(nf90_strerror(status))
          stop "Stopped by netCDF"
      endif  
    endsubroutine handle_err
    
    !! add FillValue for output
    !subroutine addFillValue(stat)
    !  type(stat_field), intent(inout) :: stat
    !  
    !  ! low left corner
    !  stat%u(ips:its-1,jps:jts-1,:) = FillValue
    !  stat%v(ips:its-1,jps:jts-1,:) = FillValue
    !  
    !  ! low right corner
    !  stat%u(ite+1:ipe,jps:jts-1,:) = FillValue
    !  stat%v(ite+1:ipe,jps:jts-1,:) = FillValue
    !  
    !  ! up left corner
    !  stat%u(ips:its-1,jte+1:jpe,:) = FillValue
    !  stat%v(ips:its-1,jte+1:jpe,:) = FillValue
    !  
    !  ! up right corner
    !  stat%u(ite+1:ipe,jte+1:jpe,:) = FillValue
    !  stat%v(ite+1:ipe,jte+1:jpe,:) = FillValue
    !end subroutine addFillValue
end module io_mod
    