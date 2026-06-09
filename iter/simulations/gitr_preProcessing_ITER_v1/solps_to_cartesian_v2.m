%% SOLPS-ITER → Cartesian grid → NetCDF export (multi-species D + Ne)
% Converts SOLPS cell data to uniform (R,Z) and writes multi-species NetCDF.

clear; clc; close all;

%% === INPUTS ===
matFile = 'solps_iter.mat';
outnc   = fullfile(pwd, 'profiles.nc');   % output path
if exist(outnc,'file'), delete(outnc); end

fprintf('Writing multiSpecies NetCDF to %s\n', outnc);

%% === LOAD DATA ===
S = load(matFile);
Geo   = S.Geo;
State = S.State;

r = double(Geo.pr);   % (4 x Nc)
z = double(Geo.pz);
Nc = size(r,2);

ne_cv = double(State.ne(:));
na_cv = double(State.na);
te_cv = double(State.te(:))./1.602e-19;
ti_cv = double(State.ti(:))./1.602e-19;
ua_cv = double(State.ua);

zn = double(State.zn(:));
am = double(State.am(:));
ns = size(na_cv,2);

fprintf('Loaded %d cells and %d species.\n', Nc, ns);

%% === DEFINE CARTESIAN GRID ===
rmin = min(r(:)); rmax = max(r(:));
zmin = min(z(:)); zmax = max(z(:));
rgrid = linspace(rmin, rmax, 400);
zgrid = linspace(zmin, zmax, 800);
[X,Z] = meshgrid(rgrid, zgrid);

fill_by_cells = @(vals) paintCells(vals, r, z, X, Z);

%% === MAP FIELDS ===
disp('Interpolating fields to (R,Z) grid ...');
ne_q = fill_by_cells(ne_cv);
Te_q = fill_by_cells(te_cv);
Ti_q = fill_by_cells(ti_cv);

% Each species density and parallel flow
ni_q = cell(1, ns);
Ua_q = cell(1, ns);
for s = 1:ns
    ni_q{s} = fill_by_cells(na_cv(:,s));
    Ua_q{s} = fill_by_cells(ua_cv(:,s));
end

%% === STACK ARRAYS ===
nR = numel(rgrid);
nZ = numel(zgrid);

ni_all = zeros(nZ, nR, ns);
ua_all = zeros(nZ, nR, ns);
Ti_all = zeros(nZ, nR, ns);

for s = 1:ns
    ni_all(:,:,s) = ni_q{s};
    ua_all(:,:,s) = Ua_q{s};
    Ti_all(:,:,s) = Ti_q;  % same Ti for all species here
end

%% === WRITE MULTISPECIES NETCDF ===
disp('Writing multiSpecies NetCDF...');

ncid = netcdf.create(outnc, bitor(netcdf.getConstant('NETCDF4'), netcdf.getConstant('CLOBBER')));

% dimensions
dimR = netcdf.defDim(ncid, 'nX', nR);
dimZ = netcdf.defDim(ncid, 'nZ', nZ);
dimS = netcdf.defDim(ncid, 'species', ns);

% variables
vid_x   = netcdf.defVar(ncid,'x','double',dimR);
vid_z   = netcdf.defVar(ncid,'z','double',dimZ);
vid_Z   = netcdf.defVar(ncid,'atomic_number','double',dimS);
vid_q   = netcdf.defVar(ncid,'charge_number','double',dimS);
vid_ne  = netcdf.defVar(ncid,'ne','double',[dimR dimZ]);
vid_te  = netcdf.defVar(ncid,'te','double',[dimR dimZ]);
vid_ti  = netcdf.defVar(ncid,'ti','double',[dimR dimZ]);
vid_niA = netcdf.defVar(ncid,'ni_all','double',[dimR dimZ dimS]);
vid_uA  = netcdf.defVar(ncid,'u_all','double',[dimR dimZ dimS]);

netcdf.endDef(ncid);

% metadata
netcdf.putVar(ncid, vid_x, rgrid);
netcdf.putVar(ncid, vid_z, zgrid);
netcdf.putVar(ncid, vid_Z, zn);
netcdf.putVar(ncid, vid_q, zn);

% single-fluid fields (transpose -> [nR x nZ])
netcdf.putVar(ncid, vid_ne, ne_q');
netcdf.putVar(ncid, vid_te, Te_q');
netcdf.putVar(ncid, vid_ti, Ti_q');

% multi-species 3D fields ([nZ x nR x ns] -> [nR x nZ x ns])
netcdf.putVar(ncid, vid_niA, permute(ni_all, [2 1 3]));
netcdf.putVar(ncid, vid_uA,  permute(ua_all, [2 1 3]));

netcdf.close(ncid);
fprintf('✅ Wrote %s (multiSpecies)\n', outnc);

%% === QUICK PLOTS ===
figure; imagesc(rgrid,zgrid,Te_q); set(gca,'YDir','normal','ColorScale','log');
colorbar; title('n_e [m^{-3}]'); xlabel('R [m]'); ylabel('Z [m]'); axis equal tight;

k = 2;
figure; imagesc(rgrid,zgrid,Ua_q{k}); set(gca,'YDir','normal','ColorScale','log');
colorbar; title(sprintf('Species %d density [m^{-3}]', k)); xlabel('R [m]'); ylabel('Z [m]');
axis equal tight;

%% === Helper ===
function grid_vals = paintCells(cv_vals, r4xN, z4xN, X, Z)
grid_vals = zeros(size(X));
N = size(r4xN,2);
for i = 1:N
    [in,on] = inpolygon(X,Z,r4xN(:,i),z4xN(:,i));
    if any(in(:)) || any(on(:))
        grid_vals(in|on) = cv_vals(i);
    end
end
end