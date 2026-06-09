close all;
clc;
clear all;
%% Read Simulation Profiles
disp('Reading simulation profiles')
x=ncread('profilesProtoMPEX.nc','x');
z=ncread('profilesProtoMPEX.nc','z');

ne=1.2.*ncread('profilesProtoMPEX.nc','ne');
figure; imagesc(z,x,ne);
set(gca,'YDir','normal')
set(gca, 'FontSize', 24)
set(gca,'FontName','times','fontSize',24);
ylabel('$r$ [m]','interpreter','Latex','fontSize',24);
xlabel('$z$ [m]','interpreter','latex','fontSize',24);
title('$n_e$', 'interpreter','latex','fontSize',26)
colorbar;
pbaspect([2 1 1])
xlim([0.5 4.2])



te=ncread('profilesProtoMPEX.nc','te');
figure; imagesc(z,x,te);
set(gca, 'FontSize', 24)
set(gca,'YDir','normal')
set(gca,'FontName','times','fontSize',24);
ylabel('$r$ [m]','interpreter','Latex','fontSize',24);
xlabel('$z$ [m]','interpreter','latex','fontSize',24);
title('$T_e$', 'interpreter','latex','fontSize',26)
colorbar;
pbaspect([2 1 1])
xlim([0.5 4.2])

% vx=ncread('../TeX1/input/profilesProtoMPEX.nc','vx');
% vx=ncread('../TeX1/input/profilesProtoMPEX.nc','vy');
vz=ncread('profilesProtoMPEX.nc','vz');
figure;imagesc(z,x,vz);
set(gca, 'FontSize', 24)
set(gca,'YDir','normal')
set(gca,'FontName','times','fontSize',24);
ylabel('$r$ [m]','interpreter','Latex','fontSize',24);
xlabel('$z$ [m]','interpreter','latex','fontSize',24);
title('$U_\parallel$', 'interpreter','latex','fontSize',26)
colorbar;
% figure; plot(vz(1,:))
pbaspect([2 1 1])
xlim([0.5 4.2])

figure; plot(z,1.5.*te(100,:),'LineWidth',2);
set(gca, 'FontSize', 24)
xlabel('$z$ [m]','interpreter','Latex','fontSize',24);
ylabel('$T_e [eV]$','interpreter','Latex','fontSize',24);
% title('$T_e$', 'interpreter','latex','fontSize',26)
figure; plot(x,te(2,2500))
xlabel('$r$ [m]','interpreter','Latex','fontSize',24);
ylabel('$T_e [eV]$','interpreter','Latex','fontSize',24);
% title('Input Radial Te')
pbaspect([2 1 1])
xlim([0.5 4.2])

figure; plot(z,1.2.*ne(2,:), 'LineWidth',2);
set(gca, 'FontSize', 24)
xlabel('$z$ [m]','interpreter','Latex','fontSize',24);
ylabel('$n_e [m^{-3}]$','interpreter','Latex','fontSize',24);
xlim([0.5 4.3])
% title('Input Axial ne')
pbaspect([2 1 1])
figure; plot(x,ne(:,2),'LineWidth',2);
set(gca, 'FontSize', 24)
xlabel('$r$ [m]','interpreter','Latex','fontSize',24);
ylabel('$n_e [m^{-3}]$','interpreter','Latex','fontSize',24);
% title('Input Radial ne')
pbaspect([2 1 1])
xlim([0.5 4.2])

figure; plot(z,vz(1,:),'LineWidth',2);
set(gca, 'FontSize', 24)
xlabel('$z$ [m]','interpreter','Latex','fontSize',24);
ylabel('$U_\parallel [m/s]$','interpreter','Latex','fontSize',24);
% title('Input Axial Vz')
pbaspect([2 1 1])
xlim([0.5 4.2])
figure; plot(x,vz(:,2500),'LineWidth',2);
set(gca, 'FontSize', 24)
xlabel('$r$ [m]','interpreter','Latex','fontSize',24);
ylabel('$U_\parallel[m/s]$','interpreter','Latex','fontSize',24);
% title('Input Radial Vz')
pbaspect([2 1 1])
xlim([0.5 4.2])






