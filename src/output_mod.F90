module output_mod
  use netcdf
  use constants_mod
  use parameters_mod
  use mesh_mod
  use stat_mod
  use projection_mod
  implicit none
    
    character(13) :: ncFile = 'mcv_output.nc'
    
    contains
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
end module output_mod
    