% clear all;close all; clc;

ME = 9.10938356e-31;
MI = 1.6737236e-27;
Q = 1.60217662e-19;
EPS0 = 8.854187e-12;

amu = 20.17; % Neon

centroid1=centroid;
data=readmatrix('plasma_profiles_highdens_equilBaseline2010.txt');
rC=sqrt(centroid1(:,1).^2+centroid1(:,2).^2);
appFit_x0 = 4.50050923; appFit_z0 = 0.49749458;
z=centroid1(:,3);

r = sqrt((rC-appFit_x0).^2+(z-appFit_z0).^2);
a= r - 3.7728 + 8.2712;
ne0=data(:,3);
te0=data(:,4);
ti0=data(:,5);
R=data(:,1);
ne=interpn(R,ne0,a);
ne(isnan(ne))=0;
te=interpn(R,te0,a);
te(isnan(te))=0;
ti=interpn(R,ti0,a);
ti(isnan(ti))=0;


vp = 0.25*sqrt(1.6E-19*te./(2*MI));
flux_surf = ne.*vp;
figure; imagesc(z,r,ne)

writematrix(ne, 'ne_surf.csv')
writematrix(te, 'te_surf.csv')
writematrix(ti, 'ti_surf.csv')
writematrix(flux_surf, 'flux_surf.csv')
