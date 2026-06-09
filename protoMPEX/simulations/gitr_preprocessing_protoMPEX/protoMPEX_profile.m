% 
% close all;
% clear all;
% clc
file = 'profilesProtoMPEX_base.nc';
% % file1='/Users/78k/Library/CloudStorage/OneDrive-OakRidgeNationalLaboratory/ORNL-ATUL-MBP/myRepos/GITR_processing/postProcessing/protoMPEX/parametricScan/no_diffision_flag2/densityScan/te8to8ne1e18to1e19/input/profilesProtoMPEX.nc';
x0 = ncread(file,'x');
z0 = ncread(file,'z');
ni0 = ncread(file,'ni');
ne0 = ncread(file,'ne');
ti0 = ncread(file,'ti');
te0 = ncread(file,'te');
gradTi0=-ncread(file,'gradTi');
vr0 = ncread(file,'vr');
vt0 = ncread(file,'vt');
vz0 = ncread(file,'vz');

% br0 = ncread(file,'br');
% bt0 = ncread(file,'bt');
% bz0 = ncread(file,'bz');
% % vz1= ncread(file1,'vz');

% figure;imagesc(x0(1:25:end-1),z0(1:25:end-1),bz0(1:25:end-1,1:25:end-1))

% Physical constants:
% =========================================================================
e_c = 1.6020e-19;
k_B = 1.3806e-23;
m_p = 1.6726e-27;
m_e = 9.1094e-31;
mu0 = 4*pi*1e-7;
c   = 299792458;
E_0 = m_p*c^2;

fileB = '../Bfield_create/bfield_protoMPEX.nc';
% file1='/Users/78k/Library/CloudStorage/OneDrive-OakRidgeNationalLaboratory/ORNL-ATUL-MBP/myRepos/GITR_processing/postProcessing/protoMPEX/parametricScan/no_diffision_flag2/densityScan/te8to8ne1e18to1e19/input/profiles_protoMPEX_SOLPS.nc';
xB = ncread(fileB,'x');
zB = ncread(fileB,'z');

br0 = ncread(fileB,'br');
bt0 = ncread(fileB,'bt');
bz0 = ncread(fileB,'bz');


x0=linspace(0,0.1,4000);
z0=linspace(0.5,4.14, 5000);

% Calculate Psi from Bz
% ----------------------
dx=xB(2)-xB(1);

for ii=1:numel(xB)
    for jj= 1:numel(zB)
        psi(ii,jj)  = 2*pi.*sum((bz0 (1:ii,jj).*xB(1:ii))).*dx;
    end
end
psi=psi';
figure;imagesc(xB(1:end-1),zB(1:end-1),psi(1:end-1,1:end-1))


%%
% Normalize Psi with the value of Psi at helicon location
% -------------------------------------------------------
[m,n]=size(psi);
% x=1:m;y=1:n;
% psi_new=(interp2(psi,0.0575/(max(x0)-min(x0))*m,1.745/(max(z0)-min(z0))*n));
psi_new=(interp2(psi,1.745/(max(zB)-min(zB))*n,0.0575/(max(xB)-min(xB))*m))


%% Nornalized Psi
% --------------
psi=psi./psi_new;

figure;imagesc(xB(1:end-1),zB(1:end-1),psi(1:end-1,1:end-1));

[rr1,zz1]=meshgrid(x0,z0);

psi_interp=interpn(zB,xB,psi,zz1,rr1);

br_interp=interpn(xB,zB,br0,rr1,zz1)';
bz_interp=interpn(xB,zB,bz0,rr1,zz1)';
bt_interp=interpn(xB,zB,bt0,rr1,zz1)';
figure;imagesc(x0(1:end-1),z0(1:end-1),psi_interp(1:end-1,1:end-1)')
psi=psi_interp;

bz_helicon=(interp2(bz_interp,1.745/(max(zB)-min(zB))*n,0.06/(max(xB)-min(xB))*m))

%% Construct a temperature profile using Psi_N for ProtoMPEX
% ==========================================================
% % Test
% te_min=2;
% te_max=2;
% 
% % Helicon
% % te_min=2;
% % te_max=8;
% % 
% ECH

% for TG
te_min=8;
te_max=2;

% % for no TG
% te_min=2;
% te_max=15;

te= abs(((psi<1).*(te_max-te_min).*(1-psi.^2).^1.75) + te_min);
te(find(psi>1))=te_max; % only for TG
figure;imagesc(z0(1:end-1),x0(1:end-1),te(1:end-1,1:end-1)');
pbaspect([3 1 1])
set(gca,'YDir','normal')
set(gca,'FontName','times','fontSize',18);
ylabel('$r$ [m]','interpreter','Latex','fontSize',18);
xlabel('$z$ [m]','interpreter','latex','fontSize',18);
title('Input Te [eV]')
colorbar;
 % figure;imagesc(x0(1:25:end-1),z0(1:25:end-1),te0(1:25:end-1,1:25:end-1)');


%{ 
%% Construct density profile using Psi_N for ProtoMPEX

% Helicon
ne_min=1.0E17;
ne_max=1.0E19;
ne= (psi<1).*(ne_max-ne_min).*(1-psi.^0.5).^1.75 + ne_min;
figure;imagesc(z0(1:end-1),x0(1:end-1),ne(1:end-1,1:end-1)');
pbaspect([3 1 1])
set(gca,'YDir','normal')
set(gca,'FontName','times','fontSize',18);
ylabel('$r$ [m]','interpreter','Latex','fontSize',18);
xlabel('$z$ [m]','interpreter','latex','fontSize',18);
title('Input ne [/m^3]')
colorbar;

%}

 %% Construct a Mach profile
% % =======================================
% 
MachNum=vz0./(sqrt((8/(2*pi))*k_B.*11604.*te0./(2*m_p)));
MachNum(isnan(MachNum))=0;
MachNum(1,:)=MachNum(2,:);
figure; plot(z0,MachNum(1,:))
% 
% 
% 
%% Update the netcdf file for GITR
% ================================

%%Modify variable here
% -------------------
ti0 = te';
te0 = te';
vz = MachNum.*(sqrt((8/(2*pi))*k_B.*11604.*te'./(2*m_p)));
% vz=vz0;
% ne0 = ne';
% ni0=ne0;
br0=br_interp;
bt0=bt_interp;
bz0=bz_interp;
vz0=vz;
figure;imagesc(z0(1:end-1),x0(1:end-1),vz0(1:end-1,1:end-1));

nR = length(x0);
nZ = length(z0);
ncid = netcdf.create('profilesProtoMPEX.nc','NC_WRITE');

dimR = netcdf.defDim(ncid,'nX',nR);
dimZ = netcdf.defDim(ncid,'nZ',nZ);

gridRnc = netcdf.defVar(ncid,'x','float',dimR);
gridZnc = netcdf.defVar(ncid,'z','float',dimZ);
Ne2Dnc = netcdf.defVar(ncid,'ne','float',[dimR dimZ]);
Ni2Dnc = netcdf.defVar(ncid,'ni','float',[dimR dimZ]);
Te2Dnc = netcdf.defVar(ncid,'te','float',[dimR dimZ]);
Ti2Dnc = netcdf.defVar(ncid,'ti','float',[dimR dimZ]);
gradTi2Dnc = netcdf.defVar(ncid,'gradTi','float',[dimR dimZ]);
vrnc = netcdf.defVar(ncid,'vr','float',[dimR dimZ]);
vtnc = netcdf.defVar(ncid,'vt','float',[dimR dimZ]);
vznc = netcdf.defVar(ncid,'vz','float',[dimR dimZ]);
brnc = netcdf.defVar(ncid,'br','float',[dimR dimZ]);
btnc = netcdf.defVar(ncid,'bt','float',[dimR dimZ]);
bznc = netcdf.defVar(ncid,'bz','float',[dimR dimZ]);

netcdf.endDef(ncid);
% 
netcdf.putVar(ncid,gridRnc,x0);
netcdf.putVar(ncid,gridZnc,z0);
netcdf.putVar(ncid,Ne2Dnc,ne0);
netcdf.putVar(ncid,Ni2Dnc,ni0);
netcdf.putVar(ncid,Te2Dnc,te0);
netcdf.putVar(ncid,Ti2Dnc,ti0);
netcdf.putVar(ncid,gradTi2Dnc,gradTi0);

netcdf.putVar(ncid,vrnc,vr0);
netcdf.putVar(ncid,vtnc,vt0);
netcdf.putVar(ncid,vznc,vz0);

netcdf.putVar(ncid,brnc,br0);
netcdf.putVar(ncid,btnc,bt0);
netcdf.putVar(ncid,bznc,bz0);

netcdf.close(ncid);
% ti1 = ncread('profilesHelicon_new.nc','ti');

% return

%% Read Simulation Profiles
disp('Reading simulation profiles')
x=ncread('profilesProtoMPEX.nc','x');
z=ncread('profilesProtoMPEX.nc','z');

ne=ncread('profilesProtoMPEX.nc','ne');
figure; imagesc(z,x,ne);
set(gca,'YDir','normal')
set(gca,'FontName','times','fontSize',24);
ylabel('$r$ [m]','interpreter','Latex','fontSize',24);
xlabel('$z$ [m]','interpreter','latex','fontSize',24);
cb = colorbar(); 
yl = ylabel(cb,'$n_e [m^{-3}]$','FontSize',20, 'Interpreter', 'latex');
pbaspect([2 1 1])




te=ncread('profilesProtoMPEX.nc','te');
figure; imagesc(z,x,te);
set(gca,'YDir','normal')
set(gca,'FontName','times','fontSize',24);
ylabel('$r$ [m]','interpreter','Latex','fontSize',24);
xlabel('$z$ [m]','interpreter','latex','fontSize',24);
cb = colorbar(); 
yl = ylabel(cb,'$T_e [eV]$','FontSize',20, 'Interpreter', 'latex');
pbaspect([2 1 1])


% vx=ncread('../TeX1/input/profiles_protoMPEX_SOLPS.nc','vx');
% vx=ncread('../TeX1/input/profiles_protoMPEX_SOLPS.nc','vy');
vz=ncread('profilesProtoMPEX.nc','vz');
figure;imagesc(z,x,vz);
set(gca,'YDir','normal')
set(gca,'FontName','times','fontSize',24);
ylabel('$r$ [m]','interpreter','Latex','fontSize',24);
xlabel('$z$ [m]','interpreter','latex','fontSize',24);
cb = colorbar(); 
yl = ylabel(cb,'$U_\parallel [m/s]$','FontSize',20, 'Interpreter', 'latex');
caxis([-1.5e4 1.5e4])
pbaspect([2 1 1])
figure; plot(z,vz(2,:))

bz=ncread('profilesProtoMPEX.nc','bz');
figure;imagesc(z,x,bz);
set(gca,'YDir','normal')
set(gca,'FontName','times','fontSize',18);
ylabel('$r$ [m]','interpreter','Latex','fontSize',18);
xlabel('$z$ [m]','interpreter','latex','fontSize',18);
title('Input Vz')
colorbar;
% figure; plot(z,bz(2,:))
% 
figure; plot(z,te(2,:)')
xlabel('$z$ [m]','interpreter','Latex','fontSize',18);
ylabel('$T_e [eV]$','interpreter','Latex','fontSize',18);
title('Input Axial Te')
figure; plot(x,te(:,2000))
xlabel('$r$ [m]','interpreter','Latex','fontSize',18);
ylabel('$T_e [eV]$','interpreter','Latex','fontSize',18);
title('Input Radial Te')

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
% 
figure; plot(z,ne(100,:)')
xlabel('$z$ [m]','interpreter','Latex','fontSize',18);
ylabel('$n_e [m^{-3}]$','interpreter','Latex','fontSize',18);
title('Input Axial ne')
figure; plot(x,ne(:,2000))
xlabel('$r$ [m]','interpreter','Latex','fontSize',18);
ylabel('$n_e [m^{-3}]$','interpreter','Latex','fontSize',18);
title('Input Radial ne')

figure; plot(z,bz(2,:))
xlabel('$z$ [m]','interpreter','Latex','fontSize',18);
ylabel('$n_e [m^{-3}]$','interpreter','Latex','fontSize',18);
title('Input Axial Bz')
% 
% figure; plot(z,vz(2500,:))
% xlabel('$z$ [m]','interpreter','Latex','fontSize',18);
% ylabel('$vz [m/s]$','interpreter','Latex','fontSize',18);
% title('Input Axial Vz')
% figure; plot(x,vz(:,1))
% xlabel('$r$ [m]','interpreter','Latex','fontSize',18);
% ylabel('$vz [m/s]$','interpreter','Latex','fontSize',18);
% title('Input Radial Vz')
% 
% 
% 
% 
% 
% 
