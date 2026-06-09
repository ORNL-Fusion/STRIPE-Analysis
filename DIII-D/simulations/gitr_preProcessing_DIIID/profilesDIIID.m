clear all;close all; clc;

m=readmatrix('Shot195196.csv');
idx=find(m(:,1)~=0);
m=m(idx,:);
%plot(m(:,1),m(:,2))
R1 = 1:0.01:3;
z1 = -1:0.01:1;
[R,z] = meshgrid(R1,z1);
ad = ((0.30912-0.52748)*R-(2.306-2.2309)*z+2.306*0.52748-0.30912*2.2309)/sqrt((0.30912-0.52748)^2+(2.306-2.2309)^2);
pcolor(R,z,ad)
shading interp
title('distance to antenna [m]')
colorbar;
xlabel('R')
ylabel('z')
vq=interp1(m(:,1),m(:,2),ad,'pchip',nan);
ne=vq;
ne(isnan(ne))=0;
figure;
pcolor(R,z,ne)
shading interp
title('density parametrized by distance to antenna [m^{-3}]')
colorbar;
xlabel('R')
ylabel('z')

ME = 9.10938356e-31;
MI = 1.6737236e-27;
Q = 1.60217662e-19;
EPS0 = 8.854187e-12;

amu = 20.17; % Neon


ni=ne;
Te=zeros(size(ne))+10;
Ti=Te;
vz=sqrt(1.6E-19.*Te./(2*MI));
vt=zeros(size(vz));
vr=zeros(size(vz));


% GITR input profiles

rgrid=R1;
zgrid=z1;

nR = length(rgrid);
nZ = length(zgrid);
ncid = netcdf.create('./profilesDIIID.nc','NC_WRITE');

dimR = netcdf.defDim(ncid,'nX',nR);
dimZ = netcdf.defDim(ncid,'nZ',nZ);

gridRnc = netcdf.defVar(ncid,'x','float',dimR);
gridZnc = netcdf.defVar(ncid,'z','float',dimZ);
Ne2Dnc = netcdf.defVar(ncid,'ne','float',[dimR dimZ]);
Ni2Dnc = netcdf.defVar(ncid,'ni','float',[dimR dimZ]);
Te2Dnc = netcdf.defVar(ncid,'te','float',[dimR dimZ]);
Ti2Dnc = netcdf.defVar(ncid,'ti','float',[dimR dimZ]);
% vrnc = netcdf.defVar(ncid,'vr','float',[dimR dimZ]);
% vtnc = netcdf.defVar(ncid,'vt','float',[dimR dimZ]);
% vznc = netcdf.defVar(ncid,'vz','float',[dimR dimZ]);
% brnc = netcdf.defVar(ncid,'br','float',[dimR dimZ]);
% btnc = netcdf.defVar(ncid,'bt','float',[dimR dimZ]);
% bznc = netcdf.defVar(ncid,'bz','float',[dimR dimZ]);

netcdf.endDef(ncid);
% 
netcdf.putVar(ncid,gridRnc,rgrid);
netcdf.putVar(ncid,gridZnc,zgrid);

netcdf.putVar(ncid,Ne2Dnc,ne');
netcdf.putVar(ncid,Ni2Dnc,ni');

netcdf.putVar(ncid,Te2Dnc,Te');
netcdf.putVar(ncid,Ti2Dnc,Ti');

% netcdf.putVar(ncid,vtnc,vt');
% netcdf.putVar(ncid,vrnc,vr');
% netcdf.putVar(ncid,vznc,vz');

%% Read GITR input data
close all; clear all;
% centroid=load('centroid_comsol.csv');
% rCentroid=sqrt(centroid(:,1).^2+centroid(:,2).^2);
% rCentroid(find(rCentroid<=3.02))=0;
R=ncread('profilesDIIID.nc','x');
z=ncread('profilesDIIID.nc','z');
% writematrix(R,'r.csv')
% writematrix(z,'z.csv')



ne=ncread('profilesDIIID.nc','ne');
figure; imagesc(R,z,ne');
set(gca,'YDir','normal')
set(gca,'FontName','times','fontSize',18);
ylabel('$r$ [m]','interpreter','Latex','fontSize',18);
xlabel('$z$ [m]','interpreter','latex','fontSize',18);
title('Input ne')
colorbar;


te=ncread('profilesDIIID.nc','te');
figure; imagesc(R,z,te');
set(gca,'YDir','normal')
set(gca,'FontName','times','fontSize',18);
ylabel('$r$ [m]','interpreter','Latex','fontSize',18);
xlabel('$z$ [m]','interpreter','latex','fontSize',18);
title('Input Te')
colorbar;

% vz=ncread('profilesDIIID.nc','vz');
% figure; imagesc(R,z,vz');
% set(gca,'YDir','normal')
% set(gca,'FontName','times','fontSize',18);
% ylabel('$r$ [m]','interpreter','Latex','fontSize',18);
% xlabel('$z$ [m]','interpreter','latex','fontSize',18);
% title('Input $V_z$', 'Interpreter','latex')
% colorbar;



