close all
clear all

% % Select directory:
% % =========================================================================

% myDir = 'LH2p0MW/output'; % west
% myDir = 'LH1p5MW/output'; % west
% myDir = 'ohmic/output'; % west
% myDir = 'gitr2022_04_28/output_100k'; % west
% myDir = 'gitr2022_05_19/output'; % west

% 
% % 
% % % Read output files:
% % % =========================================================================
% cd(myDir)
% 
% 
% highFlux=1;
% 
% erosion_rate=1.3322e+18; % LH2p0
% erosion_rate=7.5855e+17; % LH1p5
% erosion_rate=2.5712e+17; % Ohmic
% 
% % % Read position file
% % % =========================================================================
% disp('Reading the output files');
% 
file = 'positions.nc';
ncdisp(file)
x = ncread(file,'x');
y = ncread(file,'y');
z = ncread(file,'z');
vx = ncread(file,'vx');
vy = ncread(file,'vy');
vz = ncread(file,'vz');
hitWall =ncread(file,'hitWall');
weight =ncread(file,'weight');
charge =ncread(file,'charge');
nP1=length(x);

transitTime =ncread(file,'transitTime');
hit = find(hitWall);
notHit=find(hitWall==0);
notHit1 = find(hitWall==0 & weight==1);
figure(101)
histogram(weight(notHit))

% % Plot particle positions
% % ------------------------
figure
scatter3(x(hit),y(hit),z(hit))
title('Final W Particle Positions')
xlabel('X')
ylabel('Y')
zlabel('Z') 

hold on
scatter3(x(notHit),y(notHit),z(notHit))


r = sqrt(x.^2 + y.^2);
collected = find( r < 0.08 & z > 0 & z < 0.2 & hitWall > 0);

  %%  
file = strcat(pwd,'/history.nc');

x = ncread(file,'x');

y = ncread(file,'y');

z = ncread(file,'z');

vx = ncread(file,'vx');

vy = ncread(file,'vy');

vz = ncread(file,'vz');

charge = ncread(file,'charge');

weight = ncread(file,'weight');

sizeArray = size(x);

nP = sizeArray(2);
% erosionPP = erosion_rate./nP1;

hit = find(weight(end,:) < 1);

r = sqrt(x.^2 + y.^2);

% Plot patricle trajecory
 

hold on

for i=1:100:nP

scatter3(x(:,i),y(:,i),z(:,i))

end
disp('Reading g-eqdsk file');

return;
% 
% 
% 
% % % Read an g-eqdsk file
% filename= '/Users/78k/Library/CloudStorage/OneDrive-OakRidgeNationalLaboratory/ORNL-ATUL-MBP/myRepos/GITR_processing/preProcessing/west/west_54034_10p2s_mag.eqdsk';
% g = readg_g3d(filename);
% 
% % Plot poloidal magnetic field flux
% hold on
% contour(g.r,g.z,g.psirz', 'LineWidth', 1);
% plot(g.lim(1,:),g.lim(2,:),'LineWidth', 1);
% 
% set(gca,'FontName','times','fontSize',24);
% ylabel('$z$ [m]','interpreter','Latex','fontSize',24)
% xlabel('$r$ [m]','interpreter','latex','fontSize',24)
% 
% %% Ionized particle trajectories
% 
% if (exist('x1') == 0)
% 
% fid = fopen('/Users/78k/Library/CloudStorage/OneDrive-OakRidgeNationalLaboratory/ORNL-ATUL-MBP/myRepos/GITR_processing/preProcessing/west/gitrGeometryPointPlane3d.cfg');
% 
% tline = fgetl(fid);
% 
% tline = fgetl(fid);
% 
%     for i=1:18
% 
%         tline = fgetl(fid);
% 
%         evalc(tline);
% 
%     end
% 
%    Zsurface = Z;
% 
% end
% 
% figure;
% hold on;
% 
% for i=1:1:nP
% 
%  plot3(x(:,i),y(:,i),z(:,i),'k')
% 
% end
% set(gca,'FontName','times','fontSize',24);
% ylabel('$y$ [m]','interpreter','Latex','fontSize',24)
% xlabel('$x$ [m]','interpreter','latex','fontSize',24)
% zlabel('$z$ [m]','interpreter','latex','fontSize',24)
% 
% 
% surface = find(Zsurface);
% 
% nSurfaces = length(a);
% 
% %Section for finding a subset of planes to plot
% 
% r = sqrt(x.^2 + y.^2);
% 
% subset = surface;%find(r<0.07 & z1> 0.001 & z1 < .20);
% 
% % subset = find(r<0.049 & z1 > -0.001 & z1<0.001)
% 
% X = [transpose(x1(subset)),transpose(x2(subset)),transpose(x3(subset))];
% 
% Y = [transpose(y1(subset)),transpose(y2(subset)),transpose(y3(subset))];
% 
% Z = [transpose(z1(subset)),transpose(z2(subset)),transpose(z3(subset))];
% 
% %patch(transpose(X(surface,:)),transpose(Y(surface,:)),transpose(Z(surface,:)),impacts(surface),'FaceAlpha',.3)
% 
% patch(transpose(X),transpose(Y),transpose(Z),zeros(1,length(subset)),'FaceAlpha',.3,'EdgeAlpha', 0.3)%,impacts(surface);
% hold on;
% 
% for i=1:1:nP
% 
%  plot3(x(:,i),y(:,i),z(:,i),'k')
% 
% end
% set(gca,'FontName','times','fontSize',24);
% ylabel('$y$ [m]','interpreter','Latex','fontSize',24)
% xlabel('$x$ [m]','interpreter','latex','fontSize',24)
% zlabel('$z$ [m]','interpreter','latex','fontSize',24)
% 
% 
% 
% 
% file = '/Users/78k/Library/CloudStorage/OneDrive-OakRidgeNationalLaboratory/ORNL-ATUL-MBP/myRepos/GITR_processing/postProcessing/west/LH2p0MW/output/surface.nc';
% 
% grossDep0 = ncread(file,'grossDeposition');
% 
% grossEro0 = ncread(file,'grossErosion');
% % eroded_flux=readmatrix('erodedflux.csv');
% 
% % figure;
% %  patch(transpose(X),transpose(Y),transpose(Z),,'FaceAlpha',1,'EdgeAlpha', 1)%,impacts(surface)
% % 
% % 
% 
% ElementArea=area(surface)';
% 
% GrossErosion=(erosionPP.*grossEro0./ElementArea);
% 
% file = '/Users/78k/Library/CloudStorage/OneDrive-OakRidgeNationalLaboratory/ORNL-ATUL-MBP/myRepos/GITR_processing/postProcessing/west/LH2p0MW/output/surface.nc';
% 
% grossDep0 = ncread(file,'grossDeposition');
% 
% grossEro00 = ncread(file,'grossErosion');
% GrossErosion0=(erosionPP.*grossEro00./ElementArea);
% % eroded_flux=readmatrix('erodedflux.csv');
% netErosion=GrossErosion-GrossErosion0;
% 
% %  figure;
% %  patch(transpose(X),transpose(Y),transpose(Z),GrossErosion,'FaceAlpha',1,'EdgeAlpha', 1)%,impacts(surface)
% patch(transpose(X),transpose(Y),transpose(Z),GrossErosion,'FaceAlpha',1,'EdgeAlpha', 1)%,impacts(surface)
% 
% 
% %% Sheath profile
% file='/Users/78k/Library/CloudStorage/OneDrive-OakRidgeNationalLaboratory/ORNL-ATUL-MBP/myRepos/GITR_processing/postProcessing/west/LH2p0MW/output/boundary_values.nc';
% potentials = ncread(file,'potential');
% te=ncread(file,'te');
% ti=ncread(file,'ti');
% 
% subset = 1:length(x1);%find(r<0.07 & z1> 0.001 & z1 < .20);
% % subset = surface;%find(r<0.07 & z1> 0.001 & z1 < .20);
% 
% % subset = find(r<0.049 & z1 > -0.001 & z1<0.001)
% 
% X = [transpose(x1(subset)),transpose(x2(subset)),transpose(x3(subset))];
% 
% Y = [transpose(y1(subset)),transpose(y2(subset)),transpose(y3(subset))];
% 
% Z = [transpose(z1(subset)),transpose(z2(subset)),transpose(z3(subset))];
% 
% %patch(transpose(X(surface,:)),transpose(Y(surface,:)),transpose(Z(surface,:)),impacts(surface),'FaceAlpha',.3)
% 
% 
%  figure;
%  patch(transpose(X),transpose(Y),transpose(Z),potentials,'FaceAlpha',1,'EdgeAlpha', 1)%,impacts(surface)
%   figure;
%  patch(transpose(X),transpose(Y),transpose(Z),te,'FaceAlpha',1,'EdgeAlpha', 1)%,impacts(surface)
%   figure;
%  patch(transpose(X),transpose(Y),transpose(Z),ti,'FaceAlpha',1,'EdgeAlpha', 1)%,impacts(surface)
% 
 %% Density of charge states

 

specFile = 'spec.nc';

Chargedens = ncread(specFile,'n');

gridR = ncread(specFile,'gridR');

% gridY = ncread(specFile,'gridY');

gridZ = ncread(specFile,'gridZ');

 

% nP=10000;
% 
 slice1_n = Chargedens(:,:,1);
 slice1_c = Chargedens(:,:,3);

% slice1 = reshape(slice1,length(gridR),length(gridY),length(gridZ));
% 
% slice1 = sum(slice1,2);
% 
% slice1 = reshape(slice1,length(gridR),length(gridZ));

dV = (gridR(2) - gridR(1))*(gridZ(2) - gridZ(1));

dt=1e-8;

% slice1 = 1/nP*dt/dV*slice1;

figure; 
hold on;
% 
% 
% % Read an g-eqdsk file
% filename= '/Users/78k/Library/CloudStorage/OneDrive-OakRidgeNationalLaboratory/ORNL-ATUL-MBP/myRepos/WEST_CMod_data/Greg__54696/geqdsk_54696_36.4259_201x201';
% g = readg_g3d(filename);
% % imagesc(g.r, g.z, slice1');
% imagesc(gridR,gridZ,slice1');
% set(gca,'YDir','normal');
% 
% 
% % Plot poloidal magnetic field flux
% hold on
% contour(g.r,g.z,g.psirz', 'LineWidth', 2);
% plot(g.lim(1,:),g.lim(2,:),'r','LineWidth', 4);
% 
% set(gca,'FontName','times','fontSize',24);
% ylabel('$z$ [m]','interpreter','Latex','fontSize',24)
% xlabel('$r$ [m]','interpreter','latex','fontSize',24)
% 
% 
figure 
p1 = imagesc(gridR,gridZ,slice1_n');


% p.EdgeColor = 'none';
%%
% Energy distribution
% -------------------------------------------------------------------------
file = '/Users/78k/Library/CloudStorage/OneDrive-OakRidgeNationalLaboratory/ORNL-ATUL-MBP/myRepos/GITR_processing/postProcessing/west/LH2p0MW/output/surface.nc';
ncdisp(file)
edist = ncread(file,'surfEDist');
nE=100;
eDtotal = reshape(sum(edist,3),[90,nE]);

figure(105)
p4 = imagesc(linspace(0,90,90),linspace(0,1000,nE),eDtotal');
colorbar
% caxis([0 1000])
set(gca, 'YDir', 'normal');

numSurfaces=10540;

meanE = zeros(1,numSurfaces);

energy_iead_gitr = linspace(5,995,nE);

for i=1:numSurfaces

    iead = reshape(edist(:,:,i),[90,nE]);

    edist_i = sum(iead);

    meanE(i)= sum(energy_iead_gitr.*edist_i)/sum(edist_i);

end

 figure

 p4 = imagesc(linspace(0,90,90),linspace(0,100540,nE),iead');
 set(gca,'YDir','normal')

colorbar
% caxis([0 0.2])

return;



