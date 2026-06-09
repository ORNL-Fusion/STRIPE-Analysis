clc;
close all;
clear all;

%% Read GITR input data

% centroid=load('centroid_comsol.csv');
% rCentroid=sqrt(centroid(:,1).^2+centroid(:,2).^2);
% rCentroid(find(rCentroid<=3.02))=0;
R=ncread('profilesWEST.nc','x');
z=ncread('profilesWEST.nc','z');
% writematrix(R,'r.csv')
% writematrix(z,'z.csv')

bz=ncread('profilesWEST.nc','bz');
figure; imagesc(R,z,bz');
set(gca,'YDir','normal')
set(gca,'FontName','times','fontSize',18);
ylabel('$z$ [m]','interpreter','Latex','fontSize',18);
xlabel('$R$ [m]','interpreter','latex','fontSize',18);
title('Input Bz')
colorbar;

ne=ncread('profilesWEST.nc','ne');
figure; imagesc(R,z,ne');
set(gca,'YDir','normal')
set(gca,'FontName','times','fontSize',18);
ylabel('$r$ [m]','interpreter','Latex','fontSize',18);
xlabel('$z$ [m]','interpreter','latex','fontSize',18);
title('Input ne')
colorbar;
writematrix(ne,'ne.csv')

figure; plot((R),0.5.*ne(:,200));
% xlim([0 40])


% hold on,
%  scatter(rCentroid,centroid(:,3),'r')


ni=ncread('profilesWEST.nc','ni');
figure; imagesc(R,z,ni');
set(gca,'YDir','normal')
set(gca,'FontName','times','fontSize',18);
ylabel('$r$ [m]','interpreter','Latex','fontSize',18);
xlabel('$z$ [m]','interpreter','latex','fontSize',18);
title('Input ni')
colorbar;
% hold on,
%  scatter(rCentroid,centroid(:,3),'r')

% 
te=ncread('profilesWEST.nc','te');
figure; imagesc(R,z,te');
set(gca,'YDir','normal')
set(gca,'FontName','times','fontSize',18);
ylabel('$r$ [m]','interpreter','Latex','fontSize',18);
xlabel('$z$ [m]','interpreter','latex','fontSize',18);
title('Input Te')
colorbar;
writematrix(te,'te.csv')

ti=ncread('profilesWEST.nc','ti');
figure; imagesc(R,z,ti');
set(gca,'YDir','normal')
set(gca,'FontName','times','fontSize',18);
ylabel('$r$ [m]','interpreter','Latex','fontSize',18);
xlabel('$z$ [m]','interpreter','latex','fontSize',18);
title('Input Ti')
colorbar;

vt=ncread('profilesWEST.nc','vt');
figure; imagesc(R,z,(vt)');
set(gca,'YDir','normal')
set(gca,'FontName','times','fontSize',18);
ylabel('$r$ [m]','interpreter','Latex','fontSize',18);
xlabel('$z$ [m]','interpreter','latex','fontSize',18);
title('Input $V_{\parallel}$', 'Interpreter','latex')
colorbar;

vt=ncread('profilesWEST.nc','vt');
figure; imagesc(R,z,(vt)');
set(gca,'YDir','normal')
set(gca,'FontName','times','fontSize',18);
ylabel('$r$ [m]','interpreter','Latex','fontSize',18);
xlabel('$z$ [m]','interpreter','latex','fontSize',18);
title('Input $V_t$', 'Interpreter','latex')
colorbar;

vz=ncread('profilesWEST.nc','vz');
vz(isnan(vz))=0;
figure; imagesc(R,z,(vz)');
set(gca,'YDir','normal')
set(gca,'FontName','times','fontSize',18);
ylabel('$r$ [m]','interpreter','Latex','fontSize',18);
xlabel('$z$ [m]','interpreter','latex','fontSize',18);
title('Input $V_z$', 'Interpreter','latex')
colorbar;
vr=ncread('profilesWEST.nc','vr');
figure; imagesc(R,z,(vr)');
set(gca,'YDir','normal')
set(gca,'FontName','times','fontSize',18);
ylabel('$r$ [m]','interpreter','Latex','fontSize',18);
xlabel('$z$ [m]','interpreter','latex','fontSize',18);
title('Input $V_r$', 'Interpreter','latex')
colorbar;


