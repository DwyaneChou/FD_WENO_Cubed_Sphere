% This program is the test of Saturated Vapor Pressure
clc
clear

var_name = 'v';
% var_name = 'phit';
% var_name = 'zonal_wind';
% var_name = 'meridional_wind';
it       = 2;
iPatch   = 2;

nc_file = '..\run\mcv_output.nc';

dx         = ncreadatt(nc_file,'/','dx');
Fill_Value = ncreadatt(nc_file,var_name,'_FillValue');
its        = ncreadatt(nc_file,'/','its');
ite        = ncreadatt(nc_file,'/','ite');
jts        = ncreadatt(nc_file,'/','jts');
jte        = ncreadatt(nc_file,'/','jte');
ics        = ncreadatt(nc_file,'/','ics');
ice        = ncreadatt(nc_file,'/','ice');
jcs        = ncreadatt(nc_file,'/','jcs');
jce        = ncreadatt(nc_file,'/','jce');

ims = its - ics + 1;
ime = jce;
jms = jts - jcs + 1;
jme = jce;

lon      = ncread(nc_file,'lon');
lat      = ncread(nc_file,'lat');
areaCell = ncread(nc_file,'areaCell');
var      = ncread(nc_file,var_name,[1,1,1,it],[Inf,Inf,6,1]);

figure
var_plot = squeeze(var(:,:,iPatch))';
pic=pcolor(var_plot);
colormap(jet)
set(pic,'EdgeColor','None')
% shading interp

figure
var_internal = var_plot(ims:ime,jms:jme);
pic=pcolor(var_internal);
colormap(jet)
set(pic,'EdgeColor','None')
% shading interp