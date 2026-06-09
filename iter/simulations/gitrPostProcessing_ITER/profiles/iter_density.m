clear all;close all; clc;

centroid=readmatrix('centroid.csv')./1000;
data=readmatrix('plasma_profiles_highdens_equilBaseline2010.txt');
rC=sqrt(centroid(:,1).^2+centroid(:,2).^2);
appFit_x0 = 4.50050923; appFit_z0 = 0.49749458;
z=centroid(:,3);

r = sqrt((rC-appFit_x0).^2+(z-appFit_z0).^2);
a= r - 3.7728 + 8.2712;
ne=data(:,3);R=data(:,1);
den=interpn(R,ne,a);
den(isnan(den))=0;
writematrix(den, 'iter_den.csv')
% figure; patch(centroid(:,1),centroid(:,2),centroid(:,3),den
comsol_data=readmatrix('VDC4.txt');
voltage=comsol_data(:,4);
x_comsol=comsol_data(:,1);
y_comsol=comsol_data(:,2);
z_comsol=comsol_data(:,3);

[r1,z1]=meshgrid(a,z);

r_comsol=sqrt(x_comsol.^2+y_comsol.^2);
 


vol=interpn(r_comsol,z_comsol,  voltage',r1,z1);
vol(isnan(den))=0;
writematrix(vol, 'iter_den.csv')
