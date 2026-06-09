
fileSOLPS = 'interpolated_values.nc';
% fileSOLPS1='/UserS/78k/Library/CloudStorage/OneDrive-OakRidgeNationalLaboratory/ORNL-ATUL-MBP/myRepos/GITR_processing/postProcessing/protoMPEX/parametricScan/no_diffision_flag2/densityScan/te8to8ne1e18to1e19/input/profileSOLPSsProtoMPEX.nc';
rS = ncread(fileSOLPS,'gridr');
zS = ncread(fileSOLPS,'gridz');
niS = ncread(fileSOLPS,'ni');
neS = ncread(fileSOLPS,'ne');
tiS = ncread(fileSOLPS,'ti');
teS = ncread(fileSOLPS,'te');
vrS = ncread(fileSOLPS,'vr');
vtS = ncread(fileSOLPS,'vt');
vzS = ncread(fileSOLPS,'vz');
brS = ncread(fileSOLPS,'Br');
btS = ncread(fileSOLPS,'Bt');
bzS = ncread(fileSOLPS,'Bz');

file = 'profilesPISCESRF.nc';
% file1='/Users/78k/Library/CloudStorage/OneDrive-OakRidgeNationalLaboratory/ORNL-ATUL-MBP/myRepos/GITR_processing/postProcessing/protoMPEX/parametricScan/no_diffision_flag2/densityScan/te8to8ne1e18to1e19/input/profilesPISCESRF_SOLPS.nc';
x0 = ncread(file,'x');
z0 = ncread(file,'z');
ni0 = ncread(file,'ni');
ne0 = ncread(file,'ne');
ti0 = ncread(file,'ti');
te0 = ncread(file,'te');
vr0 = ncread(file,'vr');
vt0 = ncread(file,'vt');
vz0 = ncread(file,'vz');
br0 = ncread(file,'br');
bt0 = ncread(file,'bt');
bz0 = ncread(file,'bz');
% % vz1= ncread(file1,'vz');



% Physical constants:
% =========================================================================
e_c = 1.6020e-19;
k_B = 1.3806e-23;
m_p = 1.6726e-27;
m_e = 9.1094e-31;
mu0 = 4*pi*1e-7;
c   = 299792458;
E_0 = m_p*c^2;





%% density extrapolation
ne=neS;
gridr=rS;
gridz=zS;
ne_interp = 0.*ne;

figure
plot(gridr,ne(:,2500))
% for i=1:5000
%     inds = find(ne(:,i) > 0);
%     if length(inds) > 0
%         
% ne_interp_local = interp1(gridr(find(ne(:,i) > 0)), 0.16*ne(find(ne(:,i)>0),i),gridr,'linear','extrap');
% ne_interp_local(find(ne_interp_local<0)) = 0;
% ne_interp(:,i) = ne_interp_local;
%     end
% end

for i= 1:5000
    inds = find((ne(:,i)) > 0);
    if length(inds) > 0
        ne_interp(inds,i) = ne(inds,i);
        ne_interp(inds(end):end,i) = mean(ne(inds(end-5),i)).*(0.05.*(exp(-(gridr(inds(end):end) - gridr(inds(end)))./0.01)));
    end
end
ne_interp(isnan(ne_interp))=0;
figure
h = pcolor(gridz,gridr,ne_interp);
h.EdgeColor = 'none';

colorbar
title({'Proto-MPEX ne ProfileSOLPS'})
xlabel('r [m]')
ylabel('z [m]')

% figure; plot(gridz,ne(:,2500))
% hold on; plot(gridz,ne_interp(:,2500))

%% Flow profileSOLPS extrapolation
vr=vrS;
vt=vtS;
vz=vzS;
gridr=rS;
gridz=zS;



vr_interp = 0*vr;
vt_interp = 0.*vt;
vz_interp = 0.*vz;

% figure
% plot(gridr,ne(:,280))
for i= 1:5000
    inds = find(abs(vr(:,i)) > 0);
    if length(inds) > 0
        vr_interp(inds,i) = vr(inds,i);
        vr_interp(inds(end):end,i) = mean(vr(inds(end-5:end),i)).*exp(-(gridr(inds(end):end) - gridr(inds(end)))./0.01);
    end
end
vr_interp(isnan(vr_interp))=0;
for i= 1:5000
    inds = find(abs(vz(:,i)) > 0);
    if length(inds) > 0
        vz_interp(inds,i) = vz(inds,i);
        vz_interp(inds(end):end,i) = mean(vz(inds(end-5:end),i)).*exp(-(gridr(inds(end):end) - gridr(inds(end)))./0.01);
    end
end
vz_interp(isnan(vz_interp))=0;


figure
t1 = pcolor(gridz,gridr,vr_interp);
t1.EdgeColor = 'none';

title({'Proto-MPEX Vr ProfileSOLPS'})
xlabel('r [m]')
ylabel('z [m]')

figure
t2 = pcolor(gridz,gridr,vz_interp);
t2.EdgeColor = 'none';
colorbar
title({'Proto-MPEX Vz ProfileSOLPS'})
xlabel('r [m]')
ylabel('z [m]')

%% Temperature profileSOLPSs extrapolation
te=teS;
gridr=rS;
gridz=zS;
te_interp = 0.*te;

for i= 1:5000
    inds = find((te(:,i)) > 0);
    if length(inds) > 0
        te_interp(inds,i) = te(inds,i);
        te_interp(inds(end):end,i) = mean(te(inds(end-5),i)).*(1.0.*(exp(-(gridr(inds(end):end) - gridr(inds(end)))./0.01)));
    end
end

% figure
% plot(gridr,ne(:,280))
% for i=1:5000
%     inds = find(te(:,i) > 0);
%     if length(inds) > 0
% te_interp_local = interp1(gridr(find(te(:,i) > 0)), 1.4.*te(find(te(:,i)>0),i),gridr,'linear','extrap');
% te_interp_local(find(te_interp_local<0)) = 0;
% te_interp(:,i) = te_interp_local;
%     end
% end
te_interp(isnan(te_interp))=0;
figure
t = pcolor(gridz,gridr,te_interp);
t.EdgeColor = 'none';
colorbar
title({'Proto-MPEX Te ProfileSOLPS'})
xlabel('r [m]')
ylabel('z [m]')
%% Wall parameterS
ion_dens_wall = interp2(gridz, gridr, ne_interp,  1.1578, 0.0608) % z=1.1578+/- 0.1878; r=0.0608
ion_temp_wall = interp2(gridz, gridr, te_interp,  1.1578, 0.0608) % z=1.1578+/- 0.1878; r=0.0608



k=1.38e-23*11604;
c_bar = sqrt(8*k*ion_temp_wall/pi/4/1.66e-27);
flux = 0.25*ion_dens_wall*c_bar;



%% B-profile SOLPS extrapolation
% brS=br0;
% bzS=bz0;
% btS=bt0;




% br_interp = 0*br;
% bt_interp = 0.*bt;
% bz_interp = 0.*bz;
% 
% % figure
% % plot(gridr,ne(:,280))
% find(br(:,i) > 0);
% 
% 
% for i= 1:5000
%     inds = find(abs(br(:,i)) > 0);
%     if length(inds) > 0
%         br_interp(inds,i) = br(inds,i);
%         br_interp(inds(end):end,i) = mean(br(inds(end-5:end),i)).*exp(-(gridr(inds(end):end) - gridr(inds(end)))./0.01);
%     end
% end
% br_interp(isnan(br_interp))=0;
% 
% for i=1:5000
%     inds = find(bz(:,i) > 0);
%     if length(inds) > 0
% bz_interp_local = interp1(gridr(find(bz(:,i) > 0)), bz(find(bz(:,i)>0),i),gridr,'linear','extrap');
% bz_interp_local(find(bz_interp_local<0)) = 0;
% bz_interp(:,i) = bz_interp_local;
%     end
% end
% bz_interp(isnan(bz_interp))=0;
% figure
% t = pcolor(gridz,gridr,bz_interp);
% t.EdgeColor = 'none';
% colorbar
% title({'Proto-MPEX Bz ProfileSOLPS'})
% xlabel('r [m]')
% ylabel('z [m]')
% 
% 
% figure
% t1 = pcolor(gridz,gridr,br_interp);
% t1.EdgeColor = 'none';
% 
% title({'Proto-MPEX Br ProfileSOLPS'})
% xlabel('r [m]')
% ylabel('z [m]')
% 
% % figure
% % t2 = pcolor(gridz,gridr,vz_interp);
% % t2.EdgeColor = 'none';
% % colorbar
% % title({'Proto-MPEX Vz ProfileSOLPS'})
% % xlabel('r [m]')
% % ylabel('z [m]')



%% Updating SOLPS profileS in GITR
% ================================

%%Modify variable here
% -------------------
ti = te_interp;
te = te_interp;
vz = vz_interp;
vr = vr_interp;
ne = ne_interp;
ni = ne_interp;
bz = bz0;
br = br0;
bt = bt0;
% vz0=vz1;
x=rS;
z=zS;



nR = length(x);
nZ = length(z);
ncid = netcdf.create('profilesPISCESRF_SOLPS.nc','NC_WRITE');

dimR = netcdf.defDim(ncid,'nX',nR);
dimZ = netcdf.defDim(ncid,'nZ',nZ);

gridRnc = netcdf.defVar(ncid,'x','float',dimR);
gridZnc = netcdf.defVar(ncid,'z','float',dimZ);
Ne2Dnc = netcdf.defVar(ncid,'ne','float',[dimR dimZ]);
Ni2Dnc = netcdf.defVar(ncid,'ni','float',[dimR dimZ]);
Te2Dnc = netcdf.defVar(ncid,'te','float',[dimR dimZ]);
Ti2Dnc = netcdf.defVar(ncid,'ti','float',[dimR dimZ]);
vrnc = netcdf.defVar(ncid,'vr','float',[dimR dimZ]);
vtnc = netcdf.defVar(ncid,'vt','float',[dimR dimZ]);
vznc = netcdf.defVar(ncid,'vz','float',[dimR dimZ]);
brnc = netcdf.defVar(ncid,'br','float',[dimR dimZ]);
btnc = netcdf.defVar(ncid,'bt','float',[dimR dimZ]);
bznc = netcdf.defVar(ncid,'bz','float',[dimR dimZ]);

netcdf.endDef(ncid);
% 
netcdf.putVar(ncid,gridRnc,x);
netcdf.putVar(ncid,gridZnc,z);
netcdf.putVar(ncid,Ne2Dnc,ne);
netcdf.putVar(ncid,Ni2Dnc,ni);
netcdf.putVar(ncid,Te2Dnc,te);
netcdf.putVar(ncid,Ti2Dnc,ti);

netcdf.putVar(ncid,vrnc,vr);
netcdf.putVar(ncid,vtnc,vt);
netcdf.putVar(ncid,vznc,vz);

netcdf.putVar(ncid,brnc,br);
netcdf.putVar(ncid,btnc,bt);
netcdf.putVar(ncid,bznc,bz);

netcdf.close(ncid);
% ti1 = ncread('profilesHelicon_new.nc','ti');

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
% figure; plot(vz(1,:))

figure; plot(z,te(2500,:))
xlabel('$z$ [m]','interpreter','Latex','fontSize',18);
ylabel('$T_e [eV]$','interpreter','Latex','fontSize',18);
title('Input Axial Te')
figure; plot(x,te(:,1))
xlabel('$r$ [m]','interpreter','Latex','fontSize',18);
ylabel('$T_e [eV]$','interpreter','Latex','fontSize',18);
title('Input Radial Te')

figure; plot(z,ne(2500,:))
xlabel('$z$ [m]','interpreter','Latex','fontSize',18);
ylabel('$n_e [m^{-3}]$','interpreter','Latex','fontSize',18);
title('Input Axial ne')
figure; plot(x,ne(:,1))
xlabel('$r$ [m]','interpreter','Latex','fontSize',18);
ylabel('$n_e [m^{-3}]$','interpreter','Latex','fontSize',18);
title('Input Radial ne')

figure; plot(z,vz(1,:))
xlabel('$z$ [m]','interpreter','Latex','fontSize',18);
ylabel('$vz [m/s]$','interpreter','Latex','fontSize',18);
title('Input Axial Vz')
figure; plot(x,vz(:,1))
xlabel('$r$ [m]','interpreter','Latex','fontSize',18);
ylabel('$vz [m/s]$','interpreter','Latex','fontSize',18);
title('Input Radial Vz')






