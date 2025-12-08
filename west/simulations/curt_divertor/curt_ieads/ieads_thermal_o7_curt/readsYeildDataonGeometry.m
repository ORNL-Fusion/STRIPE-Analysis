close all;
clear all;
clc;

data=readmatrix('Targets.txt');
yields_data=readmatrix('yields.csv');
potential_data=data(:,1);
ne_data=data(:,2);
ni_data=data(:,11);
te_data=data(:,3);
v_data=data(:,5);
ero_data=[0;yields_data(:,8)].*ni_data.*abs(v_data);


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
patch(transpose(X),transpose(Y),transpose(Z),v_data,'FaceAlpha',1,'EdgeAlpha', 0.3)%,impacts(surface)
title('Flux')
xlabel('X [m]')
ylabel('Y [m]')
zlabel('Z [m]')
colorbar('eastoutside')

figure(5)
X = [transpose(x1(subset)),transpose(x2(subset)),transpose(x3(subset))];
Y = [transpose(y1(subset)),transpose(y2(subset)),transpose(y3(subset))];
Z = [transpose(z1(subset)),transpose(z2(subset)),transpose(z3(subset))];
%patch(transpose(X(surface,:)),transpose(Y(surface,:)),transpose(Z(surface,:)),impacts(surface),'FaceAlpha',.3)
patch(transpose(X),transpose(Y),transpose(Z),ero_data,'FaceAlpha',1,'EdgeAlpha', 0.3)%,impacts(surface)
title('Flux')
xlabel('X [m]')
ylabel('Y [m]')
zlabel('Z [m]')
colorbar('eastoutside')

centroid=readmatrix('centroid.csv');
rCentroid=sqrt(centroid(:,1).^2+centroid(:,3).^2);
rCentroid1=unique(rCentroid);
erodata1=unique(ero_data);

eroData2=interpn(ero_data1,rCentroid);




figure; plot(rCentroid, (ero_data))

r_div=[2.03902; 2.09927; 2.13635; 2.16879; 2.22904; 2.23831; 2.24758; 2.25685; 2.26612; 2.29393; 2.32174; 2.40053];
% flux=[ero_data(127); ero_data(179); ero_data(207); ero_data(235); ero_data(287); ero_data(295); ero_data(303);...
%     ero_data(311); ero_data(319); ero_data(342); ero_data(367); ero_data(425) ]

flux=[ero_data(425); ero_data(367); ero_data(342); ero_data(319); ero_data(311); ero_data(303); ...
    ero_data(295); ero_data(287); ero_data(367); ero_data(235); ero_data(207); ero_data(179) ]

 figure; plot(r_div,flux, 'bo')
