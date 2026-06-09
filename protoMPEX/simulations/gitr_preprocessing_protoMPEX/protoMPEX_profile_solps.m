
% close all;
% clear all;
% clc
% ECH cases
fileSOLPS = '../solps_protoMPEX/helicon/it1/interpolated_values.nc';
% fileSOLPS = '../solps_protoMPEX/ec_4kw/interpolated_values.nc';
% fileSOLPS = '../solps_protoMPEX/ec_8kw/interpolated_values.nc';
% fileSOLPS = '../solps_protoMPEX/ec_12kw/interpolated_values.nc';
% fileSOLPS = '../solps_protoMPEX/ec_30kw/interpolated_values.nc';
% fileSOLPS = '../solps_protoMPEX/ec_50kw/interpolated_values.nc';

% Helicon only Josh cases

% fileSOLPS = '../solps_protoMPEX/helicon_28950/interpolated_values_wTG_0.nc';
% fileSOLPS = '../solps_protoMPEX/helicon_28950/interpolated_values_woTG_0.nc';





% fileSOLPS1='/Users/78k/Library/CloudStorage/OneDrive-OakRidgeNationalLaboratory/ORNL-ATUL-MBP/data/GITR/GITR_processing/postProcessing/ProtoMPEX/analyticalScan/ECH/densityScan/Te2-30_ne1E18-1E18/input/profilesProtoMPEX_base.nc';
rS = ncread(fileSOLPS,'gridr');
zS = ncread(fileSOLPS,'gridz');
niS = ncread(fileSOLPS,'ni');
% neS = ncread(fileSOLPS,'ne'); 
neS = ncread(fileSOLPS,'ne'); % ne=10.*ne for case-6
tiS = ncread(fileSOLPS,'ti');
teS = ncread(fileSOLPS,'te');
gradTiS=-ncread(fileSOLPS,'gradTi');
gradTirS=ncread(fileSOLPS,'gradTir');
vrS = ncread(fileSOLPS,'vr');
vtS = ncread(fileSOLPS,'vt');
vzS = ncread(fileSOLPS,'vz');
brS = ncread(fileSOLPS,'Br');
btS = ncread(fileSOLPS,'Bt');
bzS = ncread(fileSOLPS,'Bz');
% 
figure; imagesc(zS,rS,neS)

set(gca,'YDir','normal')

figure; plot(zS, neS(100,:))




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
        ne_interp(inds(end):end,i) = mean(ne(inds(end-5),i)).*(0.042.*(exp(-(gridr(inds(end):end) - gridr(inds(end)))./0.01))); % case-1;
        % ne_interp(inds(end):end,i) = mean(ne(inds(end-5),i)).*(0.08.*(exp(-(gridr(inds(end):end) - gridr(inds(end)))./0.016))); % case-3
        % ne_interp(inds(end):end,i) = mean(ne(inds(end-5),i)).*(0.08.*(exp(-(gridr(inds(end):end) - gridr(inds(end)))./0.008))); % case-2;
        % ne_interp(inds(end):end,i) = mean(ne(inds(end-5),i)).*(0.5.*(exp(-(gridr(inds(end):end) - gridr(inds(end)))./0.0048))); % case-4;
        % ne_interp(inds(end):end,i) = mean(ne(inds(end-5),i)).*(1.1*(exp(-(gridr(inds(end):end) - gridr(inds(end)))./0.008))); % case-5;
        % ne_interp(inds(end):end,i) = mean(ne(inds(end-5),i)).*(0.3.*(exp(-(gridr(inds(end):end) - gridr(inds(end)))./0.08)));
        % ne_interp(inds(end):end,i) = mean(ne(inds(end-5),i)).*(1.0.*(exp(-(gridr(inds(end):end) - gridr(inds(end)))./0.0042))); % case-6; ne=10ne;
    end
end
ne_interp(isnan(ne_interp))=0;
% ne_interp(1,:)=ne_interp(2,:);
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
vr_interp(isnan(vr_interp))=0;
% vz_interp(1,:)=vz_interp(2,:);


figure
t1 = pcolor(gridz,gridr,vr_interp);
t1.EdgeColor = 'none';
colorbar

title({'Proto-MPEX Vr ProfileSOLPS'})
xlabel('r [m]')
ylabel('z [m]')

figure
t2 = pcolor(gridz,gridr,vz_interp);
t2.EdgeColor = 'none';
colorbar
title({'Proto-MPEX Vz ProfileSOLPS'})
xlabel('z [m]')
ylabel('r [m]')

%% Temperature profileSOLPSs extrapolation
te=teS;
gridr=rS;
gridz=zS;
te_interp = 0.*te;

for i= 1:5000
    inds = find((te(:,i)) > 0);
    if length(inds) > 0
        te_interp(inds,i) = te(inds,i);
        % te_interp(inds(end):end,i) = mean(te(inds(end-5),i)).*(1.0.*(exp(-(gridr(inds(end):end) - gridr(inds(end)))./0.01)));% w/TG
        te_interp(inds(end):end,i) = mean(te(inds(end-5),i)).*(1.0.*(exp(-(gridr(inds(end):end) - gridr(inds(end)))./0.01))); % w/oTG
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
% te_interp(1,:)=te_interp(2,:);
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



%{ %% B-profile SOLPS extrapolation
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
%}
%% Update coordinates

fileB = '../Bfield_create/bfield_protoMPEX.nc';
xB = double(ncread(fileB,'x'));
zB = double(ncread(fileB,'z'));

br0 = ncread(fileB,'br');
bt0 = ncread(fileB,'bt');
bz0 = ncread(fileB,'bz');

figure; imagesc(zB, xB,bz0)

% x0=linspace(0,0.0625,4000);
% z0=linspace(0.5, 4.2,5000);

[rr1,zz1]=meshgrid(rS,zS);



% ne_interp=interpn(rS,zS,ne_interp,rr1,zz1)';
% te_interp=interpn(rS,zS,te_interp,rr1,zz1)';
% vr_interp=interpn(rS,zS,vr_interp,rr1,zz1)';
% vz_interp=interpn(rS,zS,vz_interp,rr1,zz1)';
% vt_interp=interpn(rS,zS,vt_interp,rr1,zz1)';

br_interp=double(br0);
bz_interp=double(bz0);
bt_interp=double(bt0);

br_interp=interpn(xB,zB,br_interp,rr1,zz1)';
bz_interp=interpn(xB,zB,bz_interp,rr1,zz1)';
bt_interp=interpn(xB,zB,bt_interp,rr1,zz1)';

figure; imagesc(zS, rS,ne_interp)
figure; imagesc(zS, rS,ne_interp)
figure; imagesc(zS, rS,vz_interp)
figure; imagesc(zS, rS,te_interp)
figure; imagesc(zS, rS,bz_interp)


colorbar
title({'Proto-MPEX ne ProfileSOLPS'})
xlabel('r [m]')
ylabel('z [m]')


%% Updating SOLPS profileS in GITR
% ================================

%%Modify variable here
% -------------------
ti = te_interp;
te = te_interp;
vz = vz_interp;
vr = vr_interp;
vt = vt_interp;
ne = ne_interp;
ni = ne_interp;
% bz = bz0;
% br = br0;
% bt = bt0;
bz = bz_interp;
br = br_interp;
bt = bt_interp;
% % vz0=vz1;
x=rS;
z=zS;





nR = length(x);
nZ = length(z);
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
netcdf.putVar(ncid,gridRnc,x);
netcdf.putVar(ncid,gridZnc,z);
netcdf.putVar(ncid,Ne2Dnc,ne);
netcdf.putVar(ncid,Ni2Dnc,ni);
netcdf.putVar(ncid,Te2Dnc,te);
netcdf.putVar(ncid,Ti2Dnc,ti);
netcdf.putVar(ncid,gradTi2Dnc,gradTiS);

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


% figure; imagesc(z,x,gradTi./te);
% set(gca,'YDir','normal')
% set(gca,'FontName','times','fontSize',24);
% ylabel('$r$ [m]','interpreter','Latex','fontSize',24);
% xlabel('$z$ [m]','interpreter','latex','fontSize',24);
% title('Input gradTi')
% cb = colorbar(); 
% yl = ylabel(cb,'$\hat{B} \cdot \nabla T_i$','FontSize',20, 'Interpreter', 'latex');
% caxis([-10 10])
% pbaspect([2 1 1])
% 
% figure; plot(z,gradTi(2,:))
% xlabel('$z$ [m]','interpreter','Latex','fontSize',18);
% ylabel('$T_e [eV]$','interpreter','Latex','fontSize',18);
% title('Input Axial Te')
% figure; plot(z,gradTi(2400,:))
% xlabel('$r$ [m]','interpreter','Latex','fontSize',18);
% ylabel('$T_e [eV]$','interpreter','Latex','fontSize',18);
% title('Input Radial Te')
% % 
% 
% figure; plot(x,gradTi(:,4999)./te(:,4999))
% xlabel('$r$ [m]','interpreter','Latex','fontSize',18);
% ylabel('$T_e [eV]$','interpreter','Latex','fontSize',18);
% title('Input Radial Te')


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
figure; plot(z,te(2,:))
xlabel('$z$ [m]','interpreter','Latex','fontSize',18);
ylabel('$T_e [eV]$','interpreter','Latex','fontSize',18);
title('Input Axial Te')
figure; plot(z,te(2500,:))
xlabel('$r$ [m]','interpreter','Latex','fontSize',18);
ylabel('$T_e [eV]$','interpreter','Latex','fontSize',18);
title('Input Radial Te')
% 
figure; plot(z,ne(2,:))
xlabel('$z$ [m]','interpreter','Latex','fontSize',18);
ylabel('$n_e [m^{-3}]$','interpreter','Latex','fontSize',18);
title('Input Axial ne')
figure; plot(x,ne(:,2500))
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
