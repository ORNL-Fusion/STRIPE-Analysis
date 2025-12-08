clear all;close all; clc;

centroid=readmatrix('centroids.csv')./1000;
data=readmatrix('plasma_profiles_highdens_equilBaseline2010.txt');
rC=sqrt(centroid(:,1).^2+centroid(:,2).^2);
appFit_x0 = 4.50050923; appFit_z0 = 0.49749458;
z=centroid(:,3);

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
writematrix(ne, 'ne.csv')
writematrix(te, 'te.csv')
writematrix(ti, 'ti.csv')
% figure; patch(centroid(:,1),centroid(:,2),centroid(:,3),den)