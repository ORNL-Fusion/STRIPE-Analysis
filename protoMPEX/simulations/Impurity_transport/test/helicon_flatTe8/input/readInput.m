close all;
clc;
clear all;
%% Read Simulation Profiles
disp('Reading simulation profiles')
x=ncread('profilesProtoMPEX.nc','x');
z=ncread('profilesProtoMPEX.nc','z');

ne=ncread('profilesProtoMPEX.nc','ne');
figure; imagesc(z,x,ne);
set(gca,'YDir','normal')
set(gca,'FontName','times','fontSize',18);
ylabel('$r$ [m]','interpreter','Latex','fontSize',18);
xlabel('$z$ [m]','interpreter','latex','fontSize',18);
title('Input Density')
colorbar;



te=ncread('profilesProtoMPEX.nc','te');
figure; imagesc(z,x,te);
set(gca,'YDir','normal')
set(gca,'FontName','times','fontSize',18);
ylabel('$r$ [m]','interpreter','Latex','fontSize',18);
xlabel('$z$ [m]','interpreter','latex','fontSize',18);
title('Input Te')
colorbar;
gradTi=ncread('profilesProtoMPEX.nc','gradTi');
figure; imagesc(z,x,gradTi);
set(gca,'YDir','normal')
set(gca,'FontName','times','fontSize',24);
ylabel('$r$ [m]','interpreter','Latex','fontSize',24);
xlabel('$z$ [m]','interpreter','latex','fontSize',24);
title('Input gradTi')
cb = colorbar(); 
yl = ylabel(cb,'$\hat{B} \cdot \nabla T_i$','FontSize',20, 'Interpreter', 'latex');
caxis([-10 10])
pbaspect([2 1 1])


% vx=ncread('../TeX1/input/profilesProtoMPEX.nc','vx');
% vx=ncread('../TeX1/input/profilesProtoMPEX.nc','vy');
vz=ncread('profilesProtoMPEX.nc','vz');
figure;imagesc(z,x,vz);
set(gca,'YDir','normal')
set(gca,'FontName','times','fontSize',18);
ylabel('$r$ [m]','interpreter','Latex','fontSize',18);
xlabel('$z$ [m]','interpreter','latex','fontSize',18);
title('Input Vz')
colorbar;
% figure; plot(vz(1,:))

figure; plot(z,te(2,:))
xlabel('$z$ [m]','interpreter','Latex','fontSize',18);
ylabel('$T_e [eV]$','interpreter','Latex','fontSize',18);
title('Input Axial Te')
figure; plot(x,te(1,2500))
xlabel('$r$ [m]','interpreter','Latex','fontSize',18);
ylabel('$T_e [eV]$','interpreter','Latex','fontSize',18);
title('Input Radial Te')

figure; plot(z,ne(2,:))
xlabel('$z$ [m]','interpreter','Latex','fontSize',18);
ylabel('$n_e [m^{-3}]$','interpreter','Latex','fontSize',18);
title('Input Axial ne')
figure; plot(x,ne(1,2500))
xlabel('$r$ [m]','interpreter','Latex','fontSize',18);
ylabel('$n_e [m^{-3}]$','interpreter','Latex','fontSize',18);
title('Input Radial ne')

figure; plot(z,vz(2,:))
xlabel('$z$ [m]','interpreter','Latex','fontSize',18);
ylabel('$vz [m/s]$','interpreter','Latex','fontSize',18);
title('Input Axial Vz')
figure; plot(x,vz(:,1))
xlabel('$r$ [m]','interpreter','Latex','fontSize',18);
ylabel('$vz [m/s]$','interpreter','Latex','fontSize',18);
title('Input Radial Vz')






