close all;
clear all;
clc;

%% Read Plasma Profiles

fileName='psep1p5mw.h5';
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

r=reshape(r,4,165*146);
z=reshape(z,4,165*146);




zW=h5read(fileName,'/z_wall_points');
rW=h5read(fileName,'/r_wall_points');

figure(1); patch(r,z,t_o{8}(:),'EdgeColor','k');
hold on;
 
plot(rW,zW,'r-', 'LineWidth', 2);
set(gca, 'ColorScale', 'log', 'FontSize',18)
xlabel('R[m]', 'Interpreter','latex')
ylabel('z[m]', 'Interpreter','latex')
hold off;




xlim([1.5 3.5])
ylim([-1 1])
axis equal
colormap('jet')
colorbar
% 
% 
%% SOLEDGE2D to GITR grid

%% Plasma (electron + D ions) profile

% Querry points
% -------------
rq = 2.2;
zq = -0.5;
rgrid = linspace(1.8,3.3,1000);
zgrid = linspace(-1,1,2000);
rgrid=rgrid';
zgrid=zgrid';

[r_mesh z_mesh] = meshgrid(rgrid,zgrid);

% B-field profile

br_q = 0*r_mesh;
bz_q = 0*r_mesh;
bphi_q = 0*r_mesh;

    for i=1:length(r)
        [in,on] = inpolygon(r_mesh,z_mesh,r(:,i),z(:,i));
        if (length(find(in)) > 0 || length(find(on)) > 0)
        
        % figure; patch(r(:,i),z(:,i),data(i),'EdgeColor','k');
        % hold on;
        % scatter(rq,zq,'g')
        br_q(find(in)) = br(i);
        i
        end
    end
        
    for i=1:length(r)
        [in,on] = inpolygon(r_mesh,z_mesh,r(:,i),z(:,i));
        if (length(find(in)) > 0 || length(find(on)) > 0)
        
        bz_q(find(in)) = bz(i);
        i
        end
    end
        
    for i=1:length(r)
        [in,on] = inpolygon(r_mesh,z_mesh,r(:,i),z(:,i));
        if (length(find(in)) > 0 || length(find(on)) > 0)
        
        bphi_q(find(in)) = bphi(i);
        i
        end
    end

figure
h = imagesc(rgrid,zgrid,br_q);
set(gca,'YDir','normal')
% h.EdgeColor = 'none';
hold on
plot(rW,zW,'r-', 'LineWidth', 2);
colorbar
set(gca,'ColorScale','linear')
xlim([1.8 3.2])
ylim

figure
h = imagesc(rgrid,zgrid,bz_q);
set(gca,'YDir','normal')
% h.EdgeColor = 'none';
hold on
plot(rW,zW,'r-', 'LineWidth', 2);
colorbar
set(gca,'ColorScale','linear')
xlim([1.8 3.2])
ylim

figure
h = imagesc(rgrid,zgrid,bphi_q);
set(gca,'YDir','normal')
% h.EdgeColor = 'none';
hold on
plot(rW,zW,'r-', 'LineWidth', 2);
colorbar
set(gca,'ColorScale','linear')
xlim([1.8 3.2])
ylim

% Electron profile
% ----------------
ne_q = 0*r_mesh;
te_q = 0*r_mesh;
vez_q = 0*r_mesh;

for i=1:length(r)
[in,on] = inpolygon(r_mesh,z_mesh,r(:,i),z(:,i));
if (length(find(in)) > 0 || length(find(on)) > 0)

% figure; patch(r(:,i),z(:,i),data(i),'EdgeColor','k');
% hold on;
% scatter(rq,zq,'g')
ne_q(find(in)) = ne(i);
i
end
end

for i=1:length(r)
[in,on] = inpolygon(r_mesh,z_mesh,r(:,i),z(:,i));
if (length(find(in)) > 0 || length(find(on)) > 0)

te_q(find(in)) = te(i);
i
end
end

for i=1:length(r)
[in,on] = inpolygon(r_mesh,z_mesh,r(:,i),z(:,i));
if (length(find(in)) > 0 || length(find(on)) > 0)

vez_q(find(in)) = vez(i);
i
end
end

figure
h = imagesc(rgrid,zgrid,ne_q);
set(gca,'YDir','normal')
% h.EdgeColor = 'none';
hold on
plot(rW,zW,'r-', 'LineWidth', 2);
colorbar
set(gca,'ColorScale','log')
xlim([1.8 3.2])
ylim

figure
h = imagesc(rgrid,zgrid,te_q);
set(gca,'YDir','normal')
% h.EdgeColor = 'none';
hold on
plot(rW,zW,'r-', 'LineWidth', 2);
colorbar
set(gca,'ColorScale','log')
xlim([1.8 3.2])
ylim

figure
h = imagesc(rgrid,zgrid,vez_q);
set(gca,'YDir','normal')
% h.EdgeColor = 'none';
hold on
plot(rW,zW,'r-', 'LineWidth', 2);
colorbar
set(gca,'ColorScale','linear')
xlim([1.8 3.2])
ylim

%% D+ Ion profiles

% D+ Ion density
% ----------------
[r_mesh z_mesh] = meshgrid(rgrid,zgrid);
ni_q = 0*r_mesh;
ti_q = 0*r_mesh;
viz_q = 0*r_mesh;

    for i=1:length(r)
        [in,on] = inpolygon(r_mesh,z_mesh,r(:,i),z(:,i));
        if (length(find(in)) > 0 || length(find(on)) > 0)
        
        % figure; patch(r(:,i),z(:,i),data(i),'EdgeColor','k');
        % hold on;
        % scatter(rq,zq,'g')
        ni_q(find(in)) = ni(i);
        i
        end
    end
        
    for i=1:length(r)
        [in,on] = inpolygon(r_mesh,z_mesh,r(:,i),z(:,i));
        if (length(find(in)) > 0 || length(find(on)) > 0)
        
        ti_q(find(in)) = ti(i);
        i
        end
    end
        
    for i=1:length(r)
        [in,on] = inpolygon(r_mesh,z_mesh,r(:,i),z(:,i));
        if (length(find(in)) > 0 || length(find(on)) > 0)
        
        viz_q(find(in)) = viz(i);
        i
        end
    end

figure
h = imagesc(rgrid,zgrid,ni_q);
set(gca,'YDir','normal')
% h.EdgeColor = 'none';
hold on
plot(rW,zW,'r-', 'LineWidth', 2);
colorbar
set(gca,'ColorScale','log')
xlim([1.8 3.2])
ylim


figure
h = imagesc(rgrid,zgrid,ti_q);
set(gca,'YDir','normal')
% h.EdgeColor = 'none';
hold on
plot(rW,zW,'r-', 'LineWidth', 2);
colorbar
set(gca,'ColorScale','log')
xlim([1.8 3.2])
ylim

figure
h = imagesc(rgrid,zgrid,viz_q);
set(gca,'YDir','normal')
% h.EdgeColor = 'none';
hold on
plot(rW,zW,'r-', 'LineWidth', 2);
colorbar
set(gca,'ColorScale','linear')
xlim([1.8 3.2])
ylim


%% Oxygen profile

% O8+ density
% ----------------

numOxygenSpecies=8;

for ss=1:numOxygenSpecies
% O+1 profiles

speciesDensity = ['/o_',int2str(ss),'/dens'];
speciesParrFlow = ['/o_',int2str(ss),'/parr_flow'];
speciesTemp = ['/o_',int2str(ss),'/temp'];


n_o{ss} = h5read(fileName,speciesDensity)';
t_o{ss} = h5read(fileName,speciesTemp)';
vz_o{ss} = h5read(fileName,speciesParrFlow)';

no_q{ss} = 0*r_mesh;
to_q{ss} = 0*r_mesh;
vzo_q{ss} = 0*r_mesh;
    for i=1:length(r)
        [in,on] = inpolygon(r_mesh,z_mesh,r(:,i),z(:,i));
        if (length(find(in)) > 0 || length(find(on)) > 0)
    
        no_q{ss}(find(in)) = n_o{ss}(i);
        i
        end
    end

    for i=1:length(r)
        [in,on] = inpolygon(r_mesh,z_mesh,r(:,i),z(:,i));
        if (length(find(in)) > 0 || length(find(on)) > 0)
        
        to_q{ss}(find(in)) = t_o{ss}(i);
        i
        end
    end


    for i=1:length(r)
        [in,on] = inpolygon(r_mesh,z_mesh,r(:,i),z(:,i));
        if (length(find(in)) > 0 || length(find(on)) > 0)
        
        vzo_q{ss}(find(in)) = vz_o{ss}(i);
        i
        end
    end

end

figure
h = imagesc(rgrid,zgrid,no_q{8});
set(gca,'YDir','normal')
% h.EdgeColor = 'none';
hold on
plot(rW,zW,'r-', 'LineWidth', 2);
colorbar
set(gca,'ColorScale','log')
xlim([1.8 3.2])
ylim

figure
h = imagesc(rgrid,zgrid,to_q{8});
set(gca,'YDir','normal')
% h.EdgeColor = 'none';
hold on
plot(rW,zW,'r-', 'LineWidth', 2);
colorbar
set(gca,'ColorScale','linear')
xlim([1.8 3.2])
ylim

figure
h = imagesc(rgrid,zgrid,vzo_q{8});
set(gca,'YDir','normal')
% h.EdgeColor = 'none';
hold on
plot(rW,zW,'r-', 'LineWidth', 2);
colorbar
set(gca,'ColorScale','linear')
xlim([1.8 3.2])
ylim


disp('Saving .mat file for this interpolation')
fileName = fileName(1:strfind(fileName,'.')-1);
save(fileName)

%% Write data for GITR input

disp('loading SOLEDGE data')
load('psep1p5mw.mat'    )

% Initialize variables

x0=zeros(size(rgrid));
z0=zeros(size(zgrid));
br0=zeros(size(br_q));
bz0=zeros(size(bz_q));
bt0=zeros(size(bphi_q));
vt0 = zeros(size(vzo_q{8}));
vz0 = zeros(size(vzo_q{8}));
vr0 = zeros(size(vzo_q{8}));

ti0 = zeros(size(to_q{8}));
te0 = zeros(size(te_q));

ne0 = zeros(size(ne_q));
ni0=zeros(size(no_q{8}));

x0=double(rgrid);
z0=double(zgrid);

br0=double(br_q);
bt0=double(bphi_q);
bz0=double(bz_q);
ti0 = double(to_q{8});
te0 = double(te_q);
ne0 = double(ne_q);
ni0=double(no_q{8});
vt0 = double(vzo_q{8})./ni0;
vz0 = zeros(size(vzo_q{8}));
vr0 = zeros(size(vzo_q{8}));



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

%% Read Data
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

return 


%% Values at the wall
val_wall = interpn(zgrid,rgrid,ne_q,zW,rW)

soledge_wall = 0*rW;
for i=1:length(r)
[in,on] = inpolygon(rW,zW,r(:,i),z(:,i));
if (length(find(in)) > 0 || length(find(on)) > 0)


% figure; patch(r(:,i),z(:,i),data(i),'EdgeColor','k');
% hold on;
% scatter(rq,zq,'g')
soledge_wall(find(in)) = data(i);
i
end
end

figure
plot(val_wall)
hold on
plot(soledge_wall)