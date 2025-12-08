clc; clear all; close all;

yields_D = readmatrix("../ieads_D/yields_D.csv");
yields_C1 = readmatrix("../ieads_c1+/yields_c1.csv");
yields_C2 = readmatrix("../ieads_c2+/yields_c2.csv");
yields_C3 = readmatrix("../ieads_c3+/yields_c3.csv");
yields_C4 = readmatrix("../ieads_c4+/yields_c4.csv");
yields_C5 = readmatrix("../ieads_c5+/yields_c5.csv");
yields_C6 = readmatrix("../ieads_c6+/yields_c6.csv");

data_D = readmatrix("../ieads_D/Targets_D.txt");
data_C1 = readmatrix("../ieads_c1+/Targets_c1.txt");
data_C2 = readmatrix("../ieads_c2+/Targets_c2.txt");
data_C3 = readmatrix("../ieads_c3+/Targets_c3.txt");
data_C4 = readmatrix("../ieads_c4+/Targets_c4.txt");
data_C5 = readmatrix("../ieads_c5+/Targets_c5.txt");
data_C6 = readmatrix("../ieads_c6+/Targets_c6.txt");


yields_total = yields_D(:,1)+ yields_C1(:,1)+ yields_C2(:,2)+ yields_C3(:,3)+ ...
                yields_C4(:,4)+ yields_C5(:,5)+ yields_C6(:,6);
density_total = data_D(:,2)+ data_C1(:,4) + data_C2(:,4) + data_C3(:,4) + ...
                data_C4(:,4) + data_C5(:,4) + data_C6(:,4);

flow_total = [data_D(:,6)+ data_C1(:,6) + data_C2(:,6) + data_C3(:,6) + ...
                data_C4(:,6) + data_C5(:,6) + data_C6(:,6)]./7;

writematrix(yields_total,'yields_total_200882.txt');

writematrix(density_total,'density_total_200882.txt');

writematrix(flow_total,'flow_total_200882.txt');

data=readmatrix('Targets_c1.txt');
potential_data=data(:,1);
ne_data=data(:,2);
ni_data=data(:,4);
te_data=data(:,3);
v_data=data(:,6);
ero_data=yields_total.*density_total.*flow_total;
% ero_data=yields_total;


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
patch(transpose(X),transpose(Y),transpose(Z),density_total,'FaceAlpha',1,'EdgeAlpha', 0.3)%,impacts(surface)
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
patch(transpose(X),transpose(Y),transpose(Z),ero_data,'FaceAlpha',1,'EdgeAlpha', 0.3)%,impacts(surface)
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
title('Flow')
xlabel('X [m]')
ylabel('Y [m]')
zlabel('Z [m]')
colorbar('eastoutside')

% figure(5)
% erosion_rate=flux_data.*area;
% writematrix(erosion_rate,'erosion_rate.csv')
% imagesc(erosion_rate)
% title('Flux')
% xlabel('X [m]')
% ylabel('Y [m]')
% zlabel('Z [m]')
% colorbar('eastoutside')