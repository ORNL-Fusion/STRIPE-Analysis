
% clc;
% close all;
% clear all;

disp('>>>> Reading g-eqdsk file');


% Read an g-eqdsk file

% filename= 'g900003.00230_ITER_15MA_eqdsk16VVHR.txt';
filename= 'MOB-348s_eqdsk.txt';
g = readg_g3d(filename);

figure; plot_gfile(g)

set(gca,'FontName','times','fontSize',24);
ylabel('$z$ [m$^{-3}$]','interpreter','Latex','fontSize',24)
xlabel('$r$ [m]','interpreter','latex','fontSize',24)
%% Calculate B-field and write it in .csv and .nc format
disp('>>>> Writing the B-field profiles');

nR = 400; %for ITER
nZ = 400;
offset = 3;
r_efit=linspace(g.r(1+offset),g.r(end-offset),nR);
z_efit=linspace(g.z(1+offset),g.z(end-offset),nZ);
[r2D,z2D] = meshgrid(r_efit,z_efit);
[Bout,ierr] = bfield_geq_bicub(g,r2D(:),z2D(:));

[psiN,psi] = calc_psiN(g,r2D(:),z2D(:),[]);
psi=reshape(psi,[nZ,nR]);
psiN=reshape(psiN,[nZ,nR]);

Br = reshape(Bout.br,[nZ,nR]);
Bt = reshape(Bout.bphi,[nZ,nR]);
Bz = reshape(Bout.bz,[nZ,nR]);


ncid = netcdf.create(('./bfield_iter.nc'),'NC_WRITE');

dimR = netcdf.defDim(ncid,'nX',nR);

dimZ = netcdf.defDim(ncid,'nZ',nZ);

gridRnc = netcdf.defVar(ncid,'x','float',dimR);

gridZnc = netcdf.defVar(ncid,'z','float',dimZ);

brnc = netcdf.defVar(ncid,'br','float',[dimR dimZ]);
btnc = netcdf.defVar(ncid,'bt','float',[dimR dimZ]);
bznc = netcdf.defVar(ncid,'bz','float',[dimR dimZ]);

netcdf.endDef(ncid);
% 
netcdf.putVar(ncid,gridRnc,r_efit);
netcdf.putVar(ncid,gridZnc,z_efit);


netcdf.putVar(ncid,brnc,Br');
netcdf.putVar(ncid,btnc,Bt');
netcdf.putVar(ncid,bznc,Bz');

netcdf.close(ncid);

disp('>>>> Calculated the B-field profile');


