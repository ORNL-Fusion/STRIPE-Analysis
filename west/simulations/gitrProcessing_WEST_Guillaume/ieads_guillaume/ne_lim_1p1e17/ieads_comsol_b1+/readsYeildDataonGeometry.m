close all;
clear all;
clc;

data=readmatrix('Targets.txt');
yields_data=readmatrix('yields_b1+.csv');
potential_data=data(:,1);
ne_data=data(:,2);
te_data=data(:,3);
v_data=data(:,5);
ni_data=0.12.*data(:,2);
% ero_data=[0;yields_data(:,8)].*ni_data;
ero_data=abs([0;0;yields_data(:,1)].*v_data.*ne_data);




if (exist('x1') == 0)
fid = fopen(strcat('gitrGeometryPointPlane3d.cfg'));
tline = fgetl(fid);
tline = fgetl(fid);
for i=1:18
tline = fgetl(fid);
evalc(tline);
end
Zsurface = Z;
end

subset = 1:length(x1);%find(r<0.07 & z1> 0.001 & z1 < .20);
%subset = find(r<0.049 & z1 > -0.001 & z1<0.001)
figure(1)
X = [transpose(x1(subset)),transpose(x2(subset)),transpose(x3(subset))];
Y = [transpose(y1(subset)),transpose(y2(subset)),transpose(y3(subset))];
Z = [transpose(z1(subset)),transpose(z2(subset)),transpose(z3(subset))];
%patch(transpose(X(surface,:)),transpose(Y(surface,:)),transpose(Z(surface,:)),impacts(surface),'FaceAlpha',.3)
patch(transpose(X),transpose(Y),transpose(Z),potential_data,'FaceAlpha',1,'EdgeAlpha', 0.3)%,impacts(surface)
title('Yields')
colorbar('eastoutside')
xlabel('X [m]')
ylabel('Y [m]')
zlabel('Z [m]')

figure(2)
X = [transpose(x1(subset)),transpose(x2(subset)),transpose(x3(subset))];
Y = [transpose(y1(subset)),transpose(y2(subset)),transpose(y3(subset))];
Z = [transpose(z1(subset)),transpose(z2(subset)),transpose(z3(subset))];
%patch(transpose(X(surface,:)),transpose(Y(surface,:)),transpose(Z(surface,:)),impacts(surface),'FaceAlpha',.3)
patch(transpose(X),transpose(Y),transpose(Z),ne_data,'FaceAlpha',1,'EdgeAlpha', 0.3)%,impacts(surface)
title('Density')
xlabel('X [m]')
ylabel('Y [m]')
zlabel('Z [m]')
colorbar('eastoutside')

figure(3)
X = [transpose(x1(subset)),transpose(x2(subset)),transpose(x3(subset))];
Y = [transpose(y1(subset)),transpose(y2(subset)),transpose(y3(subset))];
Z = [transpose(z1(subset)),transpose(z2(subset)),transpose(z3(subset))];
%patch(transpose(X(surface,:)),transpose(Y(surface,:)),transpose(Z(surface,:)),impacts(surface),'FaceAlpha',.3)
patch(transpose(X),transpose(Y),transpose(Z),te_data,'FaceAlpha',1,'EdgeAlpha', 0.3)%,impacts(surface)
title('Temperature')
xlabel('X [m]')
ylabel('Y [m]')
zlabel('Z [m]')
colorbar('eastoutside')

figure(4)
X = [transpose(x1(subset)),transpose(x2(subset)),transpose(x3(subset))];
Y = [transpose(y1(subset)),transpose(y2(subset)),transpose(y3(subset))];
Z = [transpose(z1(subset)),transpose(z2(subset)),transpose(z3(subset))];
%patch(transpose(X(surface,:)),transpose(Y(surface,:)),transpose(Z(surface,:)),impacts(surface),'FaceAlpha',.3)
patch(transpose(X),transpose(Y),transpose(Z),ero_data,'FaceAlpha',1,'EdgeAlpha', 1)%,impacts(surface)
title('Flux')
xlabel('X [m]')
ylabel('Y [m]')
zlabel('Z [m]')
colorbar('eastoutside')
load('coolwarm.mat', 'coolwarm_rgb')
colormap(coolwarm_rgb)
set(gca, 'ColorScale', 'log');  % Apply logarithmic scale to color
clim([1e10 1e20])
axis equal


% figure(5)
% erosion_rate=flux_data.*area;
% writematrix(erosion_rate,'erosion_rate.csv')
% imagesc(erosion_rate)
% title('Flux')
% xlabel('X [m]')
% ylabel('Y [m]')
% zlabel('Z [m]')
% colorbar('eastoutside')