close all;
clear all;
clc;

data=readmatrix('Targets.txt');
% data=readmatrix('Targets_comsol_D.txt');
% data=readmatrix('Targets_comsol_o8.txt');
% data=readmatrix('Targets_comsol_o7.txt');
% data=readmatrix('Targets_comsol_o6.txt');
% data=readmatrix('Targets_comsol_o5.txt');

% data=readmatrix('Targets_thermal_o8.txt');
% data=readmatrix('Targets_thermal_o7.txt');
% data=readmatrix('Targets_thermal_o6.txt');
% data=readmatrix('Targets_thermal_o5.txt');

% data=readmatrix('Targets_thermal_o8_ohmic.txt');
% data=readmatrix('Targets_thermal_o7_ohmic.txt');
% data=readmatrix('Targets_thermal_o6_ohmic.txt');
% data=readmatrix('Targets_thermal_o5_ohmic.txt');
% 
% data=readmatrix('Targets_thermal_o8_Upar.txt');

% yields_data=readmatrix('yields_comsol_D.csv');
% yields_data=readmatrix('yields_comsol_o8.csv');
% yields_data=readmatrix('yields_comsol_o7.csv');
% yields_data=readmatrix('yields_comsol_o6.csv');
% yields_data=readmatrix('yields_comsol_o5.csv');

% yields_data=readmatrix('yields_thermal_o8.csv');
% yields_data=readmatrix('yields_thermal_o7.csv');
% yields_data=readmatrix('yields_thermal_o6.csv');
% yields_data=readmatrix('yields_thermal_o5.csv');

% yields_data=readmatrix('yields_thermal_o8_ohmic.csv');
% yields_data=readmatrix('yields_thermal_o7_ohmic.csv');
% yields_data=readmatrix('yields_thermal_o6_ohmic.csv');
% yields_data=readmatrix('yields_thermal_o5_ohmic.csv');

% yields_data=readmatrix('yields_thermal_o8_Upar.csv');



potential_data=data(:,1);
ne_data=data(:,2);
ni_data=data(:,11);
te_data=data(:,3);
ti_data=data(:,4);

k=1.38e-23*11604;
c_bar = sqrt(8*k.*ti_data/pi/4/1.66e-27);
flux = 0.25.*ni_data.*c_bar;

v_data=(data(:,5));
flux = (v_data);





% v_data(find(centroid(:,1)>=0.355))=0;
% v_data(find(centroid(:,1)<=-0.355))=0;

charge_state=8;
% ero_data=flux;
% ero_data=ni_data.*v_data;
% ero_data=0.01.*ni_data.*v_data.*[0;yields_data(:,charge_state)];

centroid=readmatrix('centroid.csv');
rCentroid=sqrt(centroid(:,1).^2+centroid(:,2).^2);

% figure; plot(rCentroid , ni_data)

if (exist('x1') == 0)
% fid = fopen(strcat('gitrGeometryPointPlane3d_comsol.cfg'));
% fid = fopen(strcat('gitrGeometryPointPlane3d_thermal.cfg'));
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
% subset = find(r<0.049 & z1 > -0.001 & z1<0.001)
figure(1)
X = [transpose(x1(subset)),transpose(x2(subset)),transpose(x3(subset))];
Y = [transpose(y1(subset)),transpose(y2(subset)),transpose(y3(subset))];
Z = [transpose(z1(subset)),transpose(z2(subset)),transpose(z3(subset))];
patch(transpose(X),transpose(Y),transpose(Z),potential_data,'FaceAlpha',1, 'EdgeAlpha', 0.3)
% patch(transpose(X),transpose(Y),transpose(Z),0.5.*[0;yields_data(:,charge_state)],'FaceAlpha',1,'EdgeAlpha', 0.3)%,impacts(surface)
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
patch(transpose(X),transpose(Y),transpose(Z),ni_data,'FaceAlpha',1,'EdgeAlpha', 0.3)%,impacts(surface)
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
title('Erosion Flux')
xlabel('X [m]')
ylabel('Y [m]')
zlabel('Z [m]')
colorbar('eastoutside')

figure(5)
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

figure(6)
X = [transpose(x1(subset)),transpose(x2(subset)),transpose(x3(subset))];
Y = [transpose(y1(subset)),transpose(y2(subset)),transpose(y3(subset))];
Z = [transpose(z1(subset)),transpose(z2(subset)),transpose(z3(subset))];
%patch(transpose(X(surface,:)),transpose(Y(surface,:)),transpose(Z(surface,:)),impacts(surface),'FaceAlpha',.3)
patch(transpose(X),transpose(Y),transpose(Z),flux,'FaceAlpha',1,'EdgeAlpha', 0.3)%,impacts(surface)
title('Upar')
xlabel('X [m]')
ylabel('Y [m]')
zlabel('Z [m]')
colorbar('eastoutside')

figure(7)
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


% int_flux=1E-4.*sum(sum(ero_data.*area))
% z_vert=[0.4; 0.3473; 0.1689; 0.0749;0.0106;-0.0851;-0.2;-0.3042;-0.4]
% 
% data_thermal=[ero_data(3863); ero_data(595); ero_data(558); ero_data(508); ero_data(4004); ero_data(926); ero_data(4044); ero_data(4091);ero_data(4141) ];
% % data_rf=[ero_data(3755); ero_data(479); ero_data(3868); ero_data(3963); ero_data(812); ero_data(3701); ero_data(4100); ero_data(4091);ero_data(4134) ];
% 
% % data=[ero_data(12741); ero_data(12773); ero_data(10864); ero_data(10888); ero_data(10942); ero_data(13590); ero_data(13573); ero_data(11072);ero_data(11096) ];
% sxbValueL=[9.58779;10.8139;12.1153;12.8246;13.4241;13.242;12.4221;11.1808;9.70449];
% % 
% % sxbValueL=[9.58779;10.8139;12.1153;12.8246;13.242;12.4221;11.1808;9.70449];
% % sxbValueR=[9.58544;10.2572;11.8979;12.9097;13.436;13.1903;12.3775;10.6288;9.92874];
% 
% % sxbValueL=[9.70449;11.1808;12.4221;13.242;12.8246;12.1153;10.8139;9.58779];
% % sxbValueR=[9.92874;10.6288;12.3775;13.436;12.9097; 11.8979;10.2572; 9.58544 ];
% % 
% % 
% % data=1E19.*[0.3026; 1.2; 1.0; 0.878; 0.4; 0.3; 0.2; 0.0 ];
% % 
% fluxL=data_thermal./(4*pi.*sxbValueL);
% % fluxR=data./sxbValueR;
% % z_vert=[0.4; 0.3;0.2;0.1;-0.1;-0.2;-0.3;-0.4]
% z_vert=[3.473; 3;1.689;0.749;-0.851;-2;-3.042;-4]
% figure; plot(z_vert,fluxR, 'ro'); hold on;
r_div=[2.03902; 2.09927; 2.13635; 2.16879; 2.22904; 2.23831; 2.24758; 2.25685; 2.26612; 2.29393; 2.32174; 2.40053];
flux=[ero_data(127); ero_data(179); ero_data(207); ero_data(235); ero_data(287); ero_data(295); ero_data(303); ero_data(311); ero_data(319); ero_data(342); ero_data(367); ero_data(425) ]
 plot(r_div,flux, 'bo')


% 

% figure(5)
% erosion_rate=flux_data.*area;
% writematrix(erosion_rate,'erosion_rate.csv')
% imagesc(erosion_rate)
% title('Flux')
% xlabel('X [m]')
% ylabel('Y [m]')
% zlabel('Z [m]')
% colorbar('eastoutside')