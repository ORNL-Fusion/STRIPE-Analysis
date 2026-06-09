close all;
clear all;
clc;

%% Read Plasma Profiles

fileName='psep1p0mw.h5';
% Grid info
% ---------
r=h5read(fileName,'/solps_like/r');
z=h5read(fileName,'/solps_like/z');

% B-field
% ---------------
br=h5read(fileName,'/bfield/b_r')';
bphi=h5read(fileName,'/bfield/b_phi')';
bz=h5read(fileName,'/bfield/b_z')';

% Electron profiles
% ----------------
ne=h5read(fileName,'/n_e/dens')';
vez=h5read(fileName,'/n_e/parr_flow')';
te=h5read(fileName,'/n_e/temp')';

% Dueterium profiles
% ------------------
ni=h5read(fileName,'/n_i/dens')';
viz=h5read(fileName,'/n_i/parr_flow')';
ti=h5read(fileName,'/n_i/temp')';

% Oxygen profiles
% ---------------

numOxygenSpecies=8;

for ss=1:numOxygenSpecies
% O+1 profiles

speciesDensity = ['/o_',int2str(ss),'/dens'];
speciesParrFlow = ['/o_',int2str(ss),'/parr_flow'];
speciesTemp = ['/o_',int2str(ss),'/temp'];


n_o{ss} = h5read(fileName,speciesDensity)';
vz_o{ss} = h5read(fileName,speciesParrFlow)';
t_o{ss} = h5read(fileName,speciesTemp)';

end

r=reshape(r,4,165*146);
z=reshape(z,4,165*146);




zW=h5read(fileName,'/z_wall_points');
rW=h5read(fileName,'/r_wall_points');

figure(1); patch(r,z,t_o{8}(:),'EdgeColor','k');
hold on;
 
plot(rW,zW,'r-', 'LineWidth', 2);
set(gca, 'ColorScale', 'log', 'FontSize',18)
xlabel('R[m]', 'Interpreter','latex')
ylabel('z[m]', 'Interpreter','latex')
hold off;




xlim([1.5 3.5])
ylim([-1 1])
axis equal
colormap('jet')
colorbar
% 
% 
%% SOLEDGE2D to GITR grid


%% Read Data
R=ncread('profilesWEST.nc','x');
z=ncread('profilesWEST.nc','z');

bz=ncread('profilesWEST.nc','bz');
figure; imagesc(R,z,bz');
set(gca,'YDir','normal')
set(gca,'FontName','times','fontSize',18);
ylabel('$z$ [m]','interpreter','Latex','fontSize',18);
xlabel('$R$ [m]','interpreter','latex','fontSize',18);
title('Input Bz')
colorbar;
hold on;
% 
plot(rW,zW,'r-', 'LineWidth', 2);
% set(gca, 'ColorScale', 'log', 'FontSize',18)
xlabel('R[m]', 'Interpreter','latex')
ylabel('z[m]', 'Interpreter','latex')
hold off;

br=ncread('profilesWEST.nc','br');
figure; imagesc(R,z,br');
set(gca,'YDir','normal')
set(gca,'FontName','times','fontSize',18);
ylabel('$z$ [m]','interpreter','Latex','fontSize',18);
xlabel('$R$ [m]','interpreter','latex','fontSize',18);
title('Input Br')
colorbar;
hold on;
% 
plot(rW,zW,'r-', 'LineWidth', 2);
% set(gca, 'ColorScale', 'log', 'FontSize',18)
xlabel('R[m]', 'Interpreter','latex')
ylabel('z[m]', 'Interpreter','latex')
hold off;

br=ncread('profilesWEST.nc','br');
figure; imagesc(R,z,br');
set(gca,'YDir','normal')
set(gca,'FontName','times','fontSize',18);
ylabel('$z$ [m]','interpreter','Latex','fontSize',18);
xlabel('$R$ [m]','interpreter','latex','fontSize',18);
title('Input Br')
colorbar;
hold on;
% 
plot(rW,zW,'r-', 'LineWidth', 2);
% set(gca, 'ColorScale', 'log', 'FontSize',18)
xlabel('R[m]', 'Interpreter','latex')
ylabel('z[m]', 'Interpreter','latex')
hold off;

bt=ncread('profilesWEST.nc','bt');
figure; imagesc(R,z,bt');
set(gca,'YDir','normal')
set(gca,'FontName','times','fontSize',18);
ylabel('$z$ [m]','interpreter','Latex','fontSize',18);
xlabel('$R$ [m]','interpreter','latex','fontSize',18);
title('Input Bt')
colorbar;
hold on;
% 
plot(rW,zW,'r-', 'LineWidth', 2);
% set(gca, 'ColorScale', 'log', 'FontSize',18)
xlabel('R[m]', 'Interpreter','latex')
ylabel('z[m]', 'Interpreter','latex')
hold off;
ne=ncread('profilesWEST.nc','ne');
figure; imagesc(R,z,ne');
set(gca,'YDir','normal')
set(gca, 'ColorScale', 'log', 'FontSize',18)
set(gca,'FontName','times','fontSize',18);
% ylabel('$r$ [m]','interpreter','Latex','fontSize',18);
% xlabel('$z$ [m]','interpreter','latex','fontSize',18);
title('Input ne')
colorbar;
hold on;
% 
plot(rW,zW,'r-', 'LineWidth', 2);
% set(gca, 'ColorScale', 'log', 'FontSize',18)
xlabel('R[m]', 'Interpreter','latex')
ylabel('z[m]', 'Interpreter','latex')

data=readmatrix("centroid.csv");
rCentroid=sqrt(data(:,1).^2+data(:,2).^2);

hold on,
 scatter(rCentroid,data(:,3),'g')
hold off



ni=ncread('profilesWEST.nc','ni');
figure; imagesc(R,z,ni');
set(gca,'YDir','normal')
set(gca,'FontName','times','fontSize',18);
ylabel('$r$ [m]','interpreter','Latex','fontSize',18);
xlabel('$z$ [m]','interpreter','latex','fontSize',18);
title('Input ni')
colorbar;

% 
te=ncread('profilesWEST.nc','te');
figure; imagesc(R,z,te');
set(gca,'YDir','normal')
set(gca,'FontName','times','fontSize',18);
% set(gca, 'ColorScale', 'log', 'FontSize',18)
ylabel('$r$ [m]','interpreter','Latex','fontSize',18);
xlabel('$z$ [m]','interpreter','latex','fontSize',18);
title('Input Te')
colorbar;

hold on;
% 
plot(rW,zW,'r-', 'LineWidth', 2);
% set(gca, 'ColorScale', 'log', 'FontSize',18)
xlabel('R[m]', 'Interpreter','latex')
ylabel('z[m]', 'Interpreter','latex')
hold off;


ti=ncread('profilesWEST.nc','ti');
figure; imagesc(R,z,ti');
set(gca,'YDir','normal')
set(gca,'FontName','times','fontSize',18);
ylabel('$r$ [m]','interpreter','Latex','fontSize',18);
xlabel('$z$ [m]','interpreter','latex','fontSize',18);
title('Input Ti')
colorbar;
hold on;
% 
plot(rW,zW,'r-', 'LineWidth', 2);
% set(gca, 'ColorScale', 'log', 'FontSize',18)
xlabel('R[m]', 'Interpreter','latex')
ylabel('z[m]', 'Interpreter','latex')
hold off;


vt=ncread('profilesWEST.nc','vt');
figure; imagesc(R,z,(vt)');
set(gca,'YDir','normal')
set(gca,'FontName','times','fontSize',18);
ylabel('$r$ [m]','interpreter','Latex','fontSize',18);
xlabel('$z$ [m]','interpreter','latex','fontSize',18);
title('Input $V_{\parallel}$', 'Interpreter','latex')
colorbar;

hold on;
% 
plot(rW,zW,'r-', 'LineWidth', 2);
% set(gca, 'ColorScale', 'log', 'FontSize',18)
xlabel('R[m]', 'Interpreter','latex')
ylabel('z[m]', 'Interpreter','latex')
hold off;


return 


%% Values at the wall
val_wall = interpn(zgrid,rgrid,ne_q,zW,rW)

soledge_wall = 0*rW;
for i=1:length(r)
[in,on] = inpolygon(rW,zW,r(:,i),z(:,i));
if (length(find(in)) > 0 || length(find(on)) > 0)


% figure; patch(r(:,i),z(:,i),data(i),'EdgeColor','k');
% hold on;
% scatter(rq,zq,'g')
soledge_wall(find(in)) = data(i);
i
end
end

figure
plot(val_wall)
hold on
plot(soledge_wall)