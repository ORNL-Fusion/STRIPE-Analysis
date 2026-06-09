%% Write data for GITR input
disp('>>>> loading SOLEDGE data')

load('psep1p5mw.mat'    )

% Initialize variables

x0=zeros(size(rgrid)      );
z0=zeros(size(zgrid)      );
br0=zeros(size(br_q)      );
bz0=zeros(size(bz_q)      );
bt0=zeros(size(bphi_q)    );
vt0 = zeros(size(vzo_q{8}));
vz0 = zeros(size(vzo_q{8}));
vr0 = zeros(size(vzo_q{8}));

ti0 = zeros(size(to_q{8}) );
te0 = zeros(size(te_q)    );

ne0 = zeros(size(ne_q)    );
ni0=zeros(size(no_q{8})   );

% Assign values to variables

x0=double(rgrid);
z0=double(zgrid);

br0=double(br_q);
bt0=double(bphi_q);
bz0=double(bz_q);

b_mag=sqrt(br0.^2+bt0.^2+bz0.^2);
% Unit vectors along B
ubr=br0./b_mag;
ubt=bt0./b_mag;
ubz=bz0./b_mag;

ti0 = double(to_q{8});
te0 = double(te_q);
ne0 = double(ne_q);
ni0=double(no_q{8});
v0 = double(vzo_q{8})./ni0;
vz0 = v0.*ubz;
vr0 = v0.*ubr;
vt0 = v0.*ubt;

vp=sqrt(vz0.^2+vt0.^2+vt0.^2);


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
netcdf.putVar(ncid,gridRnc,rgrid);
netcdf.putVar(ncid,gridZnc,zgrid);

netcdf.putVar(ncid,Ne2Dnc,ne0');
netcdf.putVar(ncid,Ni2Dnc,ni0');

netcdf.putVar(ncid,Te2Dnc,te0');
netcdf.putVar(ncid,Ti2Dnc,ti0');

netcdf.putVar(ncid,vtnc,vt0');
netcdf.putVar(ncid,vrnc,vr0');
netcdf.putVar(ncid,vznc,vz0');

netcdf.putVar(ncid,brnc,br_q');
netcdf.putVar(ncid,btnc,bphi_q');
netcdf.putVar(ncid,bznc,bz_q');
netcdf.close(ncid);

%% Read GITR input data

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

ne=ncread('profilesWEST.nc','ne');
figure; imagesc(R,z,ne');
set(gca,'YDir','normal')
set(gca,'FontName','times','fontSize',18);
ylabel('$r$ [m]','interpreter','Latex','fontSize',18);
xlabel('$z$ [m]','interpreter','latex','fontSize',18);
title('Input ne')
colorbar;


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
ylabel('$r$ [m]','interpreter','Latex','fontSize',18);
xlabel('$z$ [m]','interpreter','latex','fontSize',18);
title('Input Te')
colorbar;

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


