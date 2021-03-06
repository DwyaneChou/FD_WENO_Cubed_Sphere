% This program is the test of Saturated Vapor Pressure
clc
clear

% var_name = 'phit';
var_name = 'zonal_wind';
% var_name = 'meridional_wind';
it       = 125;

nc_file = '..\run\mcv_output.nc';

dx         = ncreadatt(nc_file,'/','dx');
Fill_Value = ncreadatt(nc_file,var_name,'_FillValue');


lon = ncread(nc_file,'lon');
lat = ncread(nc_file,'lat');
var = ncread(nc_file,var_name,[1,1,1,it],[Inf,Inf,6,1]);

var_plot = squeeze(var(:,:,5))';
pic=pcolor(var_plot);
colormap(jet)
% set(pic,'EdgeColor','None')
shading interp