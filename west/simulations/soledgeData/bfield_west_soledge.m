close all;
clear all;
clc;
load('psep1p5mw.mat');
nR=length(zgrid);
nZ=length(rgrid);
ncid = netcdf.create(('./bfield_cmod.nc'),'NC_WRITE');

dimR = netcdf.defDim(ncid,'nZ',nZ);

dimZ = netcdf.defDim(ncid,'nX',nR);

gridRnc = netcdf.defVar(ncid,'x','float',dimZ);

gridZnc = netcdf.defVar(ncid,'z','float',dimR);

brnc = netcdf.defVar(ncid,'br','float',[dimR dimZ]);
btnc = netcdf.defVar(ncid,'bt','float',[dimR dimZ]);
bznc = netcdf.defVar(ncid,'bz','float',[dimR dimZ]);

netcdf.endDef(ncid);
% 
netcdf.putVar(ncid,gridRnc,zgrid);
netcdf.putVar(ncid,gridZnc,rgrid);


netcdf.putVar(ncid,brnc,br_q');
netcdf.putVar(ncid,btnc,bphi_q');
netcdf.putVar(ncid,bznc,bz_q');

netcdf.close(ncid);

%% Read data

x=ncread('bfield_cmod.nc','x');
z=ncread('bfield_cmod.nc','z');

bz=ncread('bfield_cmod.nc','bz');
figure; imagesc(x,z,bz);
set(gca,'YDir','normal')
set(gca,'FontName','times','fontSize',18);
ylabel('$r$ [m]','interpreter','Latex','fontSize',18);
xlabel('$z$ [m]','interpreter','latex','fontSize',18);
title('Input Bz')
colorbar;

