file = 'profilesProtoMPEX.nc';
% file1='/Users/78k/Library/CloudStorage/OneDrive-OakRidgeNationalLaboratory/ORNL-ATUL-MBP/myRepos/GITR_processing/postProcessing/protoMPEX/parametricScan/no_diffision_flag2/densityScan/te8to8ne1e18to1e19/input/profilesProtoMPEX.nc';
y = ncread(file,'x');
z = ncread(file,'z');
ni = ncread(file,'ni');
ne = ncread(file,'ne');
ti = ncread(file,'ti');
te = ncread(file,'te');
vr = ncread(file,'vr');
vt = ncread(file,'vt');
vz = ncread(file,'vz');
br = ncread(file,'br');
bt = ncread(file,'bt');
bz = ncread(file,'bz');
% vz1= ncread(file1,'vz');

% Physical constants:
% =========================================================================
e_c = 1.6020e-19;
k_B = 1.3806e-23;
m_p = 1.6726e-27;
m_e = 9.1094e-31;
mu0 = 4*pi*1e-7;
c   = 299792458;
E_0 = m_p*c^2;





nY = 110;
nZ = 200;
nR = 120;
r_grid = linspace(0.0,0.08,120);
x_grid = linspace(-0.08,0.08,100);
y_grid = linspace(-0.08,0.08,110);
z_grid = linspace(-0.8+1.7462,0.8+1.7462,200);

[rr_grid,zz_grid] = meshgrid(r_grid,z_grid);


br_interp = interpn(y,z,br,rr_grid,zz_grid);
bz_interp = interpn(y,z,bz,rr_grid,zz_grid);
bt_interp = 0*br_interp;


te_interp = interpn(y,z,te,rr_grid,zz_grid);
ti_interp = interpn(y,z,ti,rr_grid,zz_grid);
ne_interp = interpn(y,z,ne,rr_grid,zz_grid);

figure(11)
h = pcolor(r_grid,z_grid,br_interp)
h.EdgeColor = 'none';
colorbar
title({'Proto-MPEX Magnetic Field (r) Profile'})
xlabel('r [m]')
ylabel('z [m]')
axis([0 max(y) min(z) max(z)])

br_interp(isnan(br_interp)) = 0;
bt_interp(isnan(bt_interp)) = 0;
bz_interp(isnan(bz_interp)) = 0;

te_interp(isnan(te_interp)) = 0;
ti_interp(isnan(ti_interp)) = 0;
ne_interp(isnan(ne_interp)) = 0;

%%

ncid = netcdf.create(['./profilesHelicon.nc'],'NC_WRITE')

dimR = netcdf.defDim(ncid,'nX',nR);
% dimY = netcdf.defDim(ncid,'nY',nY);
dimZ = netcdf.defDim(ncid,'nZ',nZ);

gridRnc = netcdf.defVar(ncid,'x','float',dimR);
% gridYnc = netcdf.defVar(ncid,'y','float',dimY);
gridZnc = netcdf.defVar(ncid,'z','float',dimZ);
Ne2Dnc = netcdf.defVar(ncid,'ne','float',[dimR dimZ]);
Ni2Dnc = netcdf.defVar(ncid,'ni','float',[dimR dimZ]);
Te2Dnc = netcdf.defVar(ncid,'te','float',[dimR dimZ]);
Ti2Dnc = netcdf.defVar(ncid,'ti','float',[dimR dimZ]);
% % vrnc = netcdf.defVar(ncid,'vr','float',[dimR dimZ]);
% % vtnc = netcdf.defVar(ncid,'vt','float',[dimR dimZ]);
% % vznc = netcdf.defVar(ncid,'vz','float',[dimR dimZ]);
brnc = netcdf.defVar(ncid,'br','float',[dimR dimZ]);
btnc = netcdf.defVar(ncid,'bt','float',[dimR dimZ]);
bznc = netcdf.defVar(ncid,'bz','float',[dimR dimZ]);
% % ernc = netcdf.defVar(ncid,'Er','float',[dimR dimZ]);
% % eznc = netcdf.defVar(ncid,'Ez','float',[dimR dimZ]);
% % etnc = netcdf.defVar(ncid,'Et','float',[dimR dimZ]);
% % gtirnc = netcdf.defVar(ncid,'gradTir','float',[dimR dimZ]);
% % gtiznc = netcdf.defVar(ncid,'gradTiz','float',[dimR dimZ]);
% % gtiync = netcdf.defVar(ncid,'gradTiy','float',[dimR dimZ]);
% % gternc = netcdf.defVar(ncid,'gradTer','float',[dimR dimZ]);
% % gteznc = netcdf.defVar(ncid,'gradTez','float',[dimR dimZ]);
% % gteync = netcdf.defVar(ncid,'gradTey','float',[dimR dimZ]);
% 
% 
% 
% %neVar = netcdf.defVar(ncid, 'Ne2', 'double',dimR);
% %teVar = netcdf.defVar(ncid, 'Te', 'double',dimR);
netcdf.endDef(ncid);
% 
netcdf.putVar(ncid,gridRnc,r_grid);
netcdf.putVar(ncid,gridZnc,z_grid);
netcdf.putVar(ncid,Ne2Dnc,ne_interp');
netcdf.putVar(ncid,Ni2Dnc,ne_interp');
netcdf.putVar(ncid,Te2Dnc,te_interp');
netcdf.putVar(ncid,Ti2Dnc,ti_interp');
% 
% % netcdf.putVar(ncid,vrnc,vrTotal');
% % netcdf.putVar(ncid,vtnc,vpTotal');
% % netcdf.putVar(ncid,vznc,vzTotal');
% 
netcdf.putVar(ncid,brnc,br_interp');
netcdf.putVar(ncid,btnc,bt_interp');
netcdf.putVar(ncid,bznc,bz_interp');
% 
% % netcdf.putVar(ncid,ernc,Epara');
% % netcdf.putVar(ncid,eznc,Eperp');
% % netcdf.putVar(ncid,etnc,(0*Epara)');
% % 
% % netcdf.putVar(ncid,gtirnc,gradTir');
% % netcdf.putVar(ncid,gtiznc,gradTiz');
% % netcdf.putVar(ncid,gtiync,gradTiy');
% % netcdf.putVar(ncid,gternc,gradTer');
% % netcdf.putVar(ncid,gteznc,gradTez');
% % netcdf.putVar(ncid,gteync,gradTey');
% 
% %netcdf.putVar(ncid, neVar, Ne2);
% %netcdf.putVar(ncid, teVar, Te);
netcdf.close(ncid);