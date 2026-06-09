% close all
% clear all;
%% Write data for GITR input
disp('>>>> loading SOLEDGE for Curt data')

fileName='plasma_background_1p5MW_o1_to_o8.h5';
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


% Initialize variables

x0=zeros(size(r));
z0=zeros(size(z));
br0=zeros(size(br));
bz0=zeros(size(bz));
bt0=zeros(size(bphi));

vte0 = zeros(size(vez));
vze0 = zeros(size(vez));
vre0 = zeros(size(vez));

vti0 = zeros(size(vez));
vzi0 = zeros(size(vez));
vri0 = zeros(size(vez));

vto0 = zeros(size(vez));
vzo0 = zeros(size(vez));
vro0 = zeros(size(vez));

ti0 = zeros(size(ti));
te0 = zeros(size(te));
to0 = zeros(size(te));

ne0 = zeros(size(ne));
ni0 = zeros(size(ni));
no0 = zeros(size(ni));

x0=double(r);
z0=double(z);

rgrid=r;
zgrid=z;

charge=1;

br0=double(br);
bt0=double(bphi);
bz0=double(bz);

te0 = double(te);
ne0 = double(ne);

ti0 = double(ti);
ni0 = double(ni);

to0 = double(t_o{charge});
no0 = double(n_o{charge});

vto0 = double(vz_o{charge})./no0;
vzo0 = zeros(size(vz_o{charge}));
vro0 = zeros(size(vz_o{charge}));

vti0 = double(viz)./ni0;
vzi0 = zeros(size(viz));
vri0 = zeros(size(viz));

vte0 = double(vez)./ne0;
vze0 = zeros(size(vez));
vre0 = zeros(size(vez));


nR = length(r);
nZ = length(z);
ncid = netcdf.create('./profilesWEST.nc','NC_WRITE');

dimR = netcdf.defDim(ncid,'nX',nR);
dimZ = netcdf.defDim(ncid,'nZ',nZ);

gridRnc = netcdf.defVar(ncid,'x','float',dimR);
gridZnc = netcdf.defVar(ncid,'z','float',dimZ);
Ne2Dnc = netcdf.defVar(ncid,'ne','float',[dimR dimZ]);
Ni2Dnc = netcdf.defVar(ncid,'ni','float',[dimR dimZ]);
No2Dnc = netcdf.defVar(ncid,'no','float',[dimR dimZ]);
Te2Dnc = netcdf.defVar(ncid,'te','float',[dimR dimZ]);
Ti2Dnc = netcdf.defVar(ncid,'ti','float',[dimR dimZ]);
To2Dnc = netcdf.defVar(ncid,'to','float',[dimR dimZ]);

vrenc = netcdf.defVar(ncid,'vre','float',[dimR dimZ]);
vtenc = netcdf.defVar(ncid,'vte','float',[dimR dimZ]);
vzenc = netcdf.defVar(ncid,'vze','float',[dimR dimZ]);

vrinc = netcdf.defVar(ncid,'vri','float',[dimR dimZ]);
vtinc = netcdf.defVar(ncid,'vti','float',[dimR dimZ]);
vzinc = netcdf.defVar(ncid,'vzi','float',[dimR dimZ]);

vronc = netcdf.defVar(ncid,'vro','float',[dimR dimZ]);
vtonc = netcdf.defVar(ncid,'vto','float',[dimR dimZ]);
vzonc = netcdf.defVar(ncid,'vzo','float',[dimR dimZ]);

brnc = netcdf.defVar(ncid,'br','float',[dimR dimZ]);
btnc = netcdf.defVar(ncid,'bt','float',[dimR dimZ]);
bznc = netcdf.defVar(ncid,'bz','float',[dimR dimZ]);

netcdf.endDef(ncid);
% 
netcdf.putVar(ncid,gridRnc,rgrid);
netcdf.putVar(ncid,gridZnc,zgrid);

netcdf.putVar(ncid,Ne2Dnc,ne0');
netcdf.putVar(ncid,Ni2Dnc,ni0');
netcdf.putVar(ncid,No2Dnc,no0');

netcdf.putVar(ncid,Te2Dnc,te0');
netcdf.putVar(ncid,Ti2Dnc,ti0');
netcdf.putVar(ncid,To2Dnc,to0');

netcdf.putVar(ncid,vtenc,vte0');
netcdf.putVar(ncid,vrenc,vre0');
netcdf.putVar(ncid,vzenc,vze0');

netcdf.putVar(ncid,vtinc,vti0');
netcdf.putVar(ncid,vrinc,vri0');
netcdf.putVar(ncid,vzinc,vzi0');

netcdf.putVar(ncid,vtonc,vto0');
netcdf.putVar(ncid,vronc,vro0');
netcdf.putVar(ncid,vzonc,vzo0');

netcdf.putVar(ncid,brnc,br');
netcdf.putVar(ncid,btnc,bphi');
netcdf.putVar(ncid,bznc,bz');
netcdf.close(ncid);

% GITR input profiles

nR = length(rgrid);
nZ = length(zgrid);
ncid = netcdf.create('./profilesWEST.nc','NC_WRITE');

dimR = netcdf.defDim(ncid,'nX',nR);
dimZ = netcdf.defDim(ncid,'nZ',nZ);

gridRnc = netcdf.defVar(ncid,'x','float',dimR);
gridZnc = netcdf.defVar(ncid,'z','float',dimZ);
Ne2Dnc = netcdf.defVar(ncid,'ne','float',[dimR dimZ]);
Ni2Dnc = netcdf.defVar(ncid,'ni','float',[dimR dimZ]);
No2Dnc = netcdf.defVar(ncid,'no','float',[dimR dimZ]);
Te2Dnc = netcdf.defVar(ncid,'te','float',[dimR dimZ]);
Ti2Dnc = netcdf.defVar(ncid,'ti','float',[dimR dimZ]);
To2Dnc = netcdf.defVar(ncid,'to','float',[dimR dimZ]);
vrenc = netcdf.defVar(ncid,'vre','float',[dimR dimZ]);
vtenc = netcdf.defVar(ncid,'vte','float',[dimR dimZ]);
vzenc = netcdf.defVar(ncid,'vze','float',[dimR dimZ]);
vrinc = netcdf.defVar(ncid,'vri','float',[dimR dimZ]);
vtinc = netcdf.defVar(ncid,'vti','float',[dimR dimZ]);
vzinc = netcdf.defVar(ncid,'vzi','float',[dimR dimZ]);
vronc = netcdf.defVar(ncid,'vro','float',[dimR dimZ]);
vtonc = netcdf.defVar(ncid,'vto','float',[dimR dimZ]);
vzonc = netcdf.defVar(ncid,'vzo','float',[dimR dimZ]);
brnc = netcdf.defVar(ncid,'br','float',[dimR dimZ]);
btnc = netcdf.defVar(ncid,'bt','float',[dimR dimZ]);
bznc = netcdf.defVar(ncid,'bz','float',[dimR dimZ]);

netcdf.endDef(ncid);
% 
netcdf.putVar(ncid,gridRnc,rgrid);
netcdf.putVar(ncid,gridZnc,zgrid);

netcdf.putVar(ncid,Ne2Dnc,ne0');
netcdf.putVar(ncid,Ni2Dnc,ni0');
netcdf.putVar(ncid,No2Dnc,no0');


netcdf.putVar(ncid,Te2Dnc,te0');
netcdf.putVar(ncid,Ti2Dnc,ti0');
netcdf.putVar(ncid,To2Dnc,to0');


netcdf.putVar(ncid,vtenc,vte0');
netcdf.putVar(ncid,vrenc,vre0');
netcdf.putVar(ncid,vzenc,vze0');

netcdf.putVar(ncid,vtinc,vti0');
netcdf.putVar(ncid,vrenc,vri0');
netcdf.putVar(ncid,vzinc,vzi0');

netcdf.putVar(ncid,vtonc,vto0');
netcdf.putVar(ncid,vronc,vro0');
netcdf.putVar(ncid,vzonc,vzo0');


netcdf.putVar(ncid,brnc,br');
netcdf.putVar(ncid,btnc,bphi');
netcdf.putVar(ncid,bznc,bz');
netcdf.close(ncid);

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

% figure; plot((R),ne(:,200));
% % xlim([0 40])

data=readmatrix("centroid.csv");
rCentroid=sqrt(data(:,1).^2+data(:,2).^2);

hold on,
 scatter(rCentroid,data(:,3),'r')


ni=ncread('profilesWEST.nc','no');
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

vt=ncread('profilesWEST.nc','vto');
figure; imagesc(R,z,(vt)');
set(gca,'YDir','normal')
set(gca,'FontName','times','fontSize',18);
ylabel('$r$ [m]','interpreter','Latex','fontSize',18);
xlabel('$z$ [m]','interpreter','latex','fontSize',18);
title('Input $V_{\parallel}$', 'Interpreter','latex')
colorbar;

vt=ncread('profilesWEST.nc','vto');
figure; imagesc(R,z,(vt)');
set(gca,'YDir','normal')
set(gca,'FontName','times','fontSize',18);
ylabel('$r$ [m]','interpreter','Latex','fontSize',18);
xlabel('$z$ [m]','interpreter','latex','fontSize',18);
title('Input $V_t$', 'Interpreter','latex')
colorbar;

vz=ncread('profilesWEST.nc','vzo');
vz(isnan(vz))=0;
figure; imagesc(R,z,(vz)');
set(gca,'YDir','normal')
set(gca,'FontName','times','fontSize',18);
ylabel('$r$ [m]','interpreter','Latex','fontSize',18);
xlabel('$z$ [m]','interpreter','latex','fontSize',18);
title('Input $V_z$', 'Interpreter','latex')
colorbar;
vr=ncread('profilesWEST.nc','vro');
figure; imagesc(R,z,(vr)');
set(gca,'YDir','normal')
set(gca,'FontName','times','fontSize',18);
ylabel('$r$ [m]','interpreter','Latex','fontSize',18);
xlabel('$z$ [m]','interpreter','latex','fontSize',18);
title('Input $V_r$', 'Interpreter','latex')
colorbar;


