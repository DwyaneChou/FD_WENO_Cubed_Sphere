% This program is the test of Saturated Vapor Pressure
% clc
% clear

% var_name = 'phit';
var_name = 'zonal_wind';
% var_name = 'meridional_wind';
it       = 2;

nc_file = '..\run\mcv_output.nc';
% nc_file = '..\run\mcv_output_RH_0p5.nc';

dx         = ncreadatt(nc_file,'/','dx');
% ids        = ncreadatt(nc_file,'/','ids');
% ide        = ncreadatt(nc_file,'/','ide');
% jds        = ncreadatt(nc_file,'/','jds');
% jde        = ncreadatt(nc_file,'/','jde');
its        = ncreadatt(nc_file,'/','its');
ite        = ncreadatt(nc_file,'/','ite');
jts        = ncreadatt(nc_file,'/','jts');
jte        = ncreadatt(nc_file,'/','jte');
ips        = ncreadatt(nc_file,'/','ips');
ipe        = ncreadatt(nc_file,'/','ipe');
jps        = ncreadatt(nc_file,'/','jps');
jpe        = ncreadatt(nc_file,'/','jpe');
ics        = ncreadatt(nc_file,'/','ics');
ice        = ncreadatt(nc_file,'/','ice');
jcs        = ncreadatt(nc_file,'/','jcs');
jce        = ncreadatt(nc_file,'/','jce');
Fill_Value = ncreadatt(nc_file,var_name,'_FillValue');

nHalo = its-ics;
Nx    = ite-its+1;
Ny    = jte-jts+1;

is    = its + nHalo;
js    = jts + nHalo;

lon = ncread(nc_file,'lon'   ,[is,js,1   ],[Nx,Ny,6  ]);
lat = ncread(nc_file,'lat'   ,[is,js,1   ],[Nx,Ny,6  ]);
var = ncread(nc_file,var_name,[is,js,1,it],[Nx,Ny,6,1]);

lon(lon<0) = 360 + lon((lon<0));

lon1d = reshape(lon,[],1);
lat1d = reshape(lat,[],1);
var1d = reshape(var,[],1);

res = dx/2;
x   = 0:res:360;
y   = -90:res:90;

[lon2d,lat2d] = meshgrid(x,y);

var_plot = griddata(lon1d,lat1d,var1d,lon2d,lat2d,'linear');

figure
pcolor(lon2d,lat2d,var_plot)
shading interp
% set(gca,'CLim',[-16,38])
colormap(jet)