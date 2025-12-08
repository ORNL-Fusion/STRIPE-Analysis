close all;
clearvars -global;
% load("extrapolated_data_196154.mat");
load("extrapolated_data_200882.mat");


%% SOLPS ITER Data on GITR Grid

%% Updating SOLPS profileS in GITR
% ================================

%%Modify variable here
% -------------------
% % 
% val_br=Br;
% val_bz=Bz;
% val_bt=Bt;
val_gradTe=val_gradTi;


vars2D = {val_ne, val_Te, val_vr, val_vz, val_vt, val_gradTi};
varNames2D = {'Electron Density (neS)', 'Electron Temperature (TeS)', ...
              'Radial Velocity (vrS)', 'Poloidal Velocity (vzS)', ...
              'Toroidal Velocity (vtS)', 'Ion Temperature Gradient (gradTiS)'};

for i = 1:length(vars2D)
    figure;
    imagesc([min(X(:)) max(X(:))], [min(Y(:)) max(Y(:))], vars2D{i});
    set(gca, 'YDir', 'normal');
    colorbar;
    set(gca, 'ColorScale', 'linear');
    title(['Extrapolated ', varNames2D{i}]);
    hold on;
    plot(g.lim(1,:), g.lim(2,:), 'r');
end

val_ne(isnan(val_ne))=0;
val_Te(isnan(val_Te))=0;
val_gradTi(isnan(val_gradTi))=0;
val_gradTe(isnan(val_gradTe))=0;

val_vr(isnan(val_vr))=0;
val_vz(isnan(val_vz))=0;
val_vt(isnan(val_vt))=0;

% val_br(isnan(val_br))=0;
% val_bz(isnan(val_bz))=0;
% val_bt(isnan(val_bt))=0;



x=X(1,:);
z=Y(:,1);

nR = length(x);
nZ = length(z);
ncid = netcdf.create('profilesDIIID_200882.nc','NC_WRITE');
% ncid = netcdf.create('profilesDIIID_196154.nc','NC_WRITE');
% ncid = netcdf.create('profilesITER.nc','NETCDF4');  % Use NetCDF-4 format

dimR = netcdf.defDim(ncid,'nX',nR);
dimZ = netcdf.defDim(ncid,'nZ',nZ);

gridRnc = netcdf.defVar(ncid,'x','float',dimR);
gridZnc = netcdf.defVar(ncid,'z','float',dimZ);
Ne2Dnc = netcdf.defVar(ncid,'ne','float',[dimR dimZ]);
Ni2Dnc = netcdf.defVar(ncid,'ni','float',[dimR dimZ]);
Te2Dnc = netcdf.defVar(ncid,'te','float',[dimR dimZ]);
Ti2Dnc = netcdf.defVar(ncid,'ti','float',[dimR dimZ]);
gradTi2Dnc = netcdf.defVar(ncid,'gradTi','float',[dimR dimZ]);
gradTe2Dnc = netcdf.defVar(ncid,'gradTe','float',[dimR dimZ]);
vrnc = netcdf.defVar(ncid,'vr','float',[dimR dimZ]);
vtnc = netcdf.defVar(ncid,'vt','float',[dimR dimZ]);
vznc = netcdf.defVar(ncid,'vz','float',[dimR dimZ]);
% brnc = netcdf.defVar(ncid,'br','float',[dimR dimZ]);
% btnc = netcdf.defVar(ncid,'bt','float',[dimR dimZ]);
% bznc = netcdf.defVar(ncid,'bz','float',[dimR dimZ]);

netcdf.endDef(ncid);
% 
netcdf.putVar(ncid,gridRnc,x);
netcdf.putVar(ncid,gridZnc,z);
netcdf.putVar(ncid,Ne2Dnc,val_ne);
netcdf.putVar(ncid,Ni2Dnc,val_ne);
netcdf.putVar(ncid,Te2Dnc,val_Te);
netcdf.putVar(ncid,Ti2Dnc,val_Te);
netcdf.putVar(ncid,gradTi2Dnc,val_gradTi);
netcdf.putVar(ncid,gradTe2Dnc,val_gradTe);

netcdf.putVar(ncid,vrnc,val_vr);
netcdf.putVar(ncid,vtnc,val_vt);
netcdf.putVar(ncid,vznc,val_vz);

% netcdf.putVar(ncid,brnc,val_br);
% netcdf.putVar(ncid,btnc,val_bt);
% netcdf.putVar(ncid,bznc,val_bz);

netcdf.close(ncid);
% ti1 = ncread('profilesHelicon_new.nc','ti');

%% Read Simulation Profiles
disp('Reading simulation profiles')
x=ncread('profilesDIIID_200882.nc','x');
z=ncread('profilesDIIID_200882.nc','z');

x=ncread('profilesDIIID_196154.nc','x');
z=ncread('profilesDIIID_196154.nc','z');

% ne=ncread('profilesDIIID_200882.nc','ne');
ne=ncread('profilesDIIID_196154.nc','ne');
figure; imagesc(z,x,ne);
set(gca,'YDir','normal')
set(gca,'FontName','times','fontSize',24);
ylabel('$r$ [m]','interpreter','Latex','fontSize',24);
xlabel('$z$ [m]','interpreter','latex','fontSize',24);
cb = colorbar(); 
yl = ylabel(cb,'$n_e [m^{-3}]$','FontSize',20, 'Interpreter', 'latex');
pbaspect([2 1 1])




% te=ncread('profilesDIIID_200882.nc','te');
te=ncread('profilesDIIID_196154.nc','te');
figure; imagesc(z,x,te);
set(gca,'YDir','normal')
set(gca,'FontName','times','fontSize',24);
ylabel('$r$ [m]','interpreter','Latex','fontSize',24);
xlabel('$z$ [m]','interpreter','latex','fontSize',24);
cb = colorbar(); 
yl = ylabel(cb,'$T_e [eV]$','FontSize',20, 'Interpreter', 'latex');
pbaspect([2 1 1])

% gradTi=ncread('profilesDIIID_200882.nc','gradTi');
gradTi=ncread('profilesDIIID_196154.nc','gradTi');
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



% vz=ncread('profilesDIIID_200882.nc','vz');
vz=ncread('profilesDIIID_196154.nc','vz');
figure;imagesc(z,x,vz);
set(gca,'YDir','normal')
set(gca,'FontName','times','fontSize',24);
ylabel('$r$ [m]','interpreter','Latex','fontSize',24);
xlabel('$z$ [m]','interpreter','latex','fontSize',24);
cb = colorbar(); 
yl = ylabel(cb,'$U_\parallel [m/s]$','FontSize',20, 'Interpreter', 'latex');
caxis([-1.5e4 1.5e4])
pbaspect([2 1 1])
figure; plot(z,vz(500,:))


% 
figure; plot(z,te(2,:))
xlabel('$z$ [m]','interpreter','Latex','fontSize',18);
ylabel('$T_e [eV]$','interpreter','Latex','fontSize',18);
title('Input Axial Te')
figure; plot(z,te(500,:))
xlabel('$r$ [m]','interpreter','Latex','fontSize',18);
ylabel('$T_e [eV]$','interpreter','Latex','fontSize',18);
title('Input Radial Te')
% 
figure; plot(z,ne(2,:))
xlabel('$z$ [m]','interpreter','Latex','fontSize',18);
ylabel('$n_e [m^{-3}]$','interpreter','Latex','fontSize',18);
title('Input Axial ne')
figure; plot(x,ne(:,500))
xlabel('$r$ [m]','interpreter','Latex','fontSize',18);
ylabel('$n_e [m^{-3}]$','interpreter','Latex','fontSize',18);
title('Input Radial ne')

