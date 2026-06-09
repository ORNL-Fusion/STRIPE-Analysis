close all;
clc;
clear all;
%% Read Simulation Profiles
disp('Reading simulation profiles')
x=ncread('profilesPISCESRF_SOLPS.nc','x');
z=ncread('profilesPISCESRF_SOLPS.nc','z');

ne=ncread('profilesPISCESRF_SOLPS.nc','ne');
figure; imagesc(z,x,ne);
set(gca,'YDir','normal')
set(gca,'FontName','times','fontSize',18);
ylabel('$r$ [m]','interpreter','Latex','fontSize',18);
xlabel('$z$ [m]','interpreter','latex','fontSize',18);
title('Input Density')
colorbar;



te=ncread('profilesPISCESRF_SOLPS.nc','te');
figure; imagesc(z,x,te);
set(gca,'YDir','normal')
set(gca,'FontName','times','fontSize',18);
ylabel('$r$ [m]','interpreter','Latex','fontSize',18);
xlabel('$z$ [m]','interpreter','latex','fontSize',18);
title('Input Te')
colorbar;

% vx=ncread('../TeX1/input/profilesPISCESRF_SOLPS_case1.nc','vx');
% vx=ncread('../TeX1/input/profilesPISCESRF_SOLPS_case1.nc','vy');
vz=ncread('profilesPISCESRF.nc','vz');
figure;imagesc(z,x,vz);
set(gca,'YDir','normal')
set(gca,'FontName','times','fontSize',18);
ylabel('$r$ [m]','interpreter','Latex','fontSize',18);
xlabel('$z$ [m]','interpreter','latex','fontSize',18);
title('Input Vz')
colorbar;
% figure; plot(vz(1,:))

figure; plot(z,te(2500,:))
xlabel('$z$ [m]','interpreter','Latex','fontSize',18);
ylabel('$T_e [eV]$','interpreter','Latex','fontSize',18);
title('Input Axial Te')
figure; plot(x,te(:,2500))
xlabel('$r$ [m]','interpreter','Latex','fontSize',18);
ylabel('$T_e [eV]$','interpreter','Latex','fontSize',18);
title('Input Radial Te')

figure; plot(z,ne(2500,:))
xlabel('$z$ [m]','interpreter','Latex','fontSize',18);
ylabel('$n_e [m^{-3}]$','interpreter','Latex','fontSize',18);
title('Input Axial ne')
figure; plot(x,ne(:,2500))
xlabel('$r$ [m]','interpreter','Latex','fontSize',18);
ylabel('$n_e [m^{-3}]$','interpreter','Latex','fontSize',18);
title('Input Radial ne')

figure; plot(z,vz(2500,:))
xlabel('$z$ [m]','interpreter','Latex','fontSize',18);
ylabel('$vz [m/s]$','interpreter','Latex','fontSize',18);
title('Input Axial Vz')
figure; plot(x,vz(:,2500))
xlabel('$r$ [m]','interpreter','Latex','fontSize',18);
ylabel('$vz [m/s]$','interpreter','Latex','fontSize',18);
title('Input Radial Vz')






