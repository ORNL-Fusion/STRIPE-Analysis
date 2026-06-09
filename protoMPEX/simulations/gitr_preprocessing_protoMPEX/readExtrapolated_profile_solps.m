close all;
clear all;
clc;

x0=readmatrix('avec.csv')';
z0=readmatrix('zvec.csv')';
bz0=readmatrix('bz.csv');
br0=readmatrix('br.csv');
bt0=zeros(size(bz0));
bz0=abs(bz0);
% % z0=z0-0.34;



% Physical constants:
% =========================================================================
e_c = 1.6020e-19;
k_B = 1.3806e-23;
m_p = 1.6726e-27;
m_e = 9.1094e-31;
mu0 = 4*pi*1e-7;
c   = 299792458;
E_0 = m_p*c^2;


% Calculate Psi from Bz
% ----------------------
dx=x0(2)-x0(1);

for ii=1:numel(x0)
    for jj= 1:numel(z0)
        psi(ii,jj)  = 2*pi.*sum((bz0 (1:ii,jj).*x0(1:ii))).*dx;
    end
end

% Normalize Psi with the value of Psi at helicon location
% -------------------------------------------------------
[m,n]=size(psi);
% x=1:m;y=1:n;
% psi_new=(interp2(psi,0.0575/(max(x0)-min(x0))*m,1.745/(max(z0)-min(z0))*n));
psi_new=(interp2(psi,1.1578/(max(z0)-min(z0))*n,0.0603/(max(x0)-min(x0))*m));


% Nornalized Psi
% --------------
psi=psi./psi_new;
% psi(psi<0)=0;



%% Read Simulation Profiles
disp('Reading simulation profiles')
x=ncread('profilesPISCESRF_SOLPS.nc','x');
z=ncread('profilesPISCESRF_SOLPS.nc','z');

bz=ncread('profilesPISCESRF_SOLPS.nc','bz');
bz0=ncread('profilesPISCESRF_SOLPS.nc','bz');
rS = ncread('interpolated_values.nc','gridr');
zS = ncread('interpolated_values.nc','gridz');
figure; imagesc(z,x,bz);
set(gca,'YDir','normal')
set(gca,'FontName','times','fontSize',18);
ylabel('$r$ [m]','interpreter','Latex','fontSize',18);
xlabel('$z$ [m]','interpreter','latex','fontSize',18);
title('Bz')
colorbar;

figure;
yyaxis left 
imagesc(z,x,bz0);
set(gca,'YDir','normal')
hold on;
yyaxis right
plot(z,bz0(2,:),'k-', 'LineWidth', 2); hold on;
plot(z,bz(2,:),'r', 'LineWidth', 4); 



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

% vx=ncread('../TeX1/input/profilesPISCESRF_SOLPS.nc','vx');
% vx=ncread('../TeX1/input/profilesPISCESRF_SOLPS.nc','vy');
vz=ncread('profilesPISCESRF_SOLPS.nc','vz');
figure;imagesc(z,x,vz);
set(gca,'YDir','normal')
set(gca,'FontName','times','fontSize',18);
ylabel('$r$ [m]','interpreter','Latex','fontSize',18);
xlabel('$z$ [m]','interpreter','latex','fontSize',18);
title('Input Vz')
colorbar;

figure;plot(z,vz(2,:));
figure;plot(z,ne(2,:));
figure;plot(z,te(2,:));

% ion_dens_wall = interp2(z, x, ne,  1.1578, 0.0608) % z=1.1578+/- 0.1878; r=0.0608
% ion_temp_wall = interp2(z, x, te,  1.1578, 0.0608) % z=1.1578+/- 0.1878; r=0.0608






