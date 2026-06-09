%% ITER Antenna Surface Plasma Profiles + B-field (with patch X,Y,Z)
% ======================================================================
% Combines:
%   - SOLPS ITER plasma profiles (profiles.nc)
%   - EFIT B-field (bfield_iter.nc)
%   - ITER antenna geometry (iter_geom.mat)
% Produces per-face plasma quantities and smooth 3-D patch plots.
% ======================================================================

clear; clc; close all;

%% -------------------- Load geometry --------------------
load('final_ITERGeom.mat','planes','centroid','norm_vec','X','Y','Z');
fprintf('Loaded ITER geometry: %d triangular faces\n', size(planes,1));

r_centroid = sqrt(centroid(:,1).^2 + centroid(:,2).^2);
phi_centroid = atan2(centroid(:,2),centroid(:,1));

%% -------------------- Load SOLPS plasma profiles --------------------
solps_file = 'profiles.nc';
fprintf('Reading SOLPS plasma profiles from %s\n', solps_file);

Rp = ncread(solps_file,'x');  Zp = ncread(solps_file,'z');
[Rp_grid, Zp_grid] = meshgrid(Rp, Zp);
ne = ncread(solps_file,'ne');  te = ncread(solps_file,'te');  ti = ncread(solps_file,'ti');
ni_all = ncread(solps_file,'ni_all');
Z_all = ncread(solps_file,'atomic_number');  q_all = ncread(solps_file,'charge_number');
[nR,nZ,nS] = size(ni_all);
fprintf('SOLPS grid: nR=%d, nZ=%d, nspecies=%d\n', nR,nZ,nS);

% Robust species mapping (handles missing metadata)
idx_D1  = find(Z_all==1  & q_all==1,1);
idx_Ne0 = find(Z_all==10 & q_all==0,1);
idx_Ne10= find(Z_all==10 & q_all==10,1);
if isempty(idx_D1)||isempty(idx_Ne0)||isempty(idx_Ne10)
    warning('Species metadata incomplete — using fallback [D1=2, Ne0=3, Ne10=13].');
    idx_D1=2; idx_Ne0=3; idx_Ne10=min(13,nS);
end
fprintf('Using species: D+ = %d, Ne0 = %d, Ne10+ = %d\n',idx_D1,idx_Ne0,idx_Ne10);

%% -------------------- Load EFIT B-field --------------------
bfile = 'bfield_iter.nc';
fprintf('Reading B-field from %s\n', bfile);

Rb = ncread(bfile,'x');  Zb = ncread(bfile,'z');
[Rb_grid,Zb_grid] = meshgrid(Rb,Zb);
Br = ncread(bfile,'br')';  Bt = ncread(bfile,'bt')';  Bz = ncread(bfile,'bz')';

%% -------------------- Interpolate profiles to centroids --------------------
interp_method = 'natural';
swap = @(A) permute(A,[2 1]);  % helper to swap R–Z

ne_surf     = interp2(Rp_grid,Zp_grid,swap(ne),   r_centroid,centroid(:,3),interp_method);
te_surf     = interp2(Rp_grid,Zp_grid,swap(te),   r_centroid,centroid(:,3),interp_method);
ti_surf     = interp2(Rp_grid,Zp_grid,swap(ti),   r_centroid,centroid(:,3),interp_method);
ni_D_surf   = interp2(Rp_grid,Zp_grid,swap(ni_all(:,:,idx_D1)),  r_centroid,centroid(:,3),interp_method);
ni_Ne0_surf = interp2(Rp_grid,Zp_grid,swap(ni_all(:,:,idx_Ne0)), r_centroid,centroid(:,3),interp_method);
ni_Ne10_surf= interp2(Rp_grid,Zp_grid,swap(ni_all(:,:,idx_Ne10)),r_centroid,centroid(:,3),interp_method);

Br_surf = interp2(Rb_grid,Zb_grid,Br,r_centroid,centroid(:,3),interp_method);
Bt_surf = interp2(Rb_grid,Zb_grid,Bt,r_centroid,centroid(:,3),interp_method);
Bz_surf = interp2(Rb_grid,Zb_grid,Bz,r_centroid,centroid(:,3),interp_method);

% Fill NaNs
vars={'ne_surf','te_surf','ti_surf','ni_D_surf','ni_Ne0_surf','ni_Ne10_surf','Br_surf','Bt_surf','Bz_surf'};
for i=1:numel(vars); eval([vars{i} '(isnan(' vars{i} '))=0;']); end

%% -------------------- B-field: toroidal → Cartesian --------------------
bx = Br_surf.*cos(phi_centroid) - Bt_surf.*sin(phi_centroid);
by = Br_surf.*sin(phi_centroid) + Bt_surf.*cos(phi_centroid);
bz = Bz_surf;
b_mag = sqrt(bx.^2 + by.^2 + bz.^2);
ubx=bx./b_mag; uby=by./b_mag; ubz=bz./b_mag;

%% -------------------- Surface normals & θ(B,n) --------------------
norm_mag = sqrt(sum(norm_vec.^2,2));
unorm = norm_vec ./ norm_mag;
theta = acos(unorm(:,1).*ubx + unorm(:,2).*uby + unorm(:,3).*ubz);
theta(isnan(theta))=0;  theta(theta>pi/2)=pi-theta(theta>pi/2);

%% -------------------- Ne10+ flux --------------------
Ne10_flux = ni_Ne10_surf .* b_mag;
writematrix(Ne10_flux,'ne10_flux_surf.csv');
fprintf('✅ Saved Ne10+ flux: ne10_flux_surf.csv\n');

%% -------------------- 3-D patch visualization (fast, per-face colors) --------------------
fprintf('Plotting ITER antenna surface fields (fast flat-shaded)...\n');

% ---- Build a single patch mesh from planes ----
nF = size(planes,1);
VA = planes(:,1:3);
VB = planes(:,4:6);
VC = planes(:,7:9);
V  = [VA; VB; VC];                              % all vertices stacked
F  = [(1:nF)' (1:nF)'+nF (1:nF)'+2*nF];         % face indices (nF×3)

% ---- Quantities to plot (one scalar per face) ----
Cset = {ne_surf, ni_D_surf, ni_Ne0_surf, ni_Ne10_surf, te_surf, ti_surf};
titles = {'Electron Density n_e [m^{-3}]', ...
          'Deuterium Ion Density n_{D^+} [m^{-3}]', ...
          'Neon Neutral Density n_{Ne^0} [m^{-3}]', ...
          'Neon Ion Density n_{Ne^{10+}} [m^{-3}]', ...
          'Electron Temperature T_e [eV]', ...
          'Ion Temperature T_i [eV]'};

for k = 1:numel(Cset)
    figure('Color','w');
    patch('Faces',F,'Vertices',V,...
          'FaceVertexCData',Cset{k},...   % one color per face
          'FaceColor','flat','EdgeColor','none');
    axis equal tight vis3d
    xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]');
    title(titles{k},'FontSize',15);
    colorbar; colormap(parula);
    view(35,25); 
end


%% -------------------- Diagnostics --------------------
fprintf('\n=== ITER Antenna Diagnostics ===\n');
fprintf('n_e: %.2e–%.2e [m^-3]\n',min(ne_surf),max(ne_surf));
fprintf('n_D^+: %.2e–%.2e [m^-3]\n',min(ni_D_surf),max(ni_D_surf));
fprintf('n_Ne^0: %.2e–%.2e [m^-3]\n',min(ni_Ne0_surf),max(ni_Ne0_surf));
fprintf('n_Ne^{10+}: %.2e–%.2e [m^-3]\n',min(ni_Ne10_surf),max(ni_Ne10_surf));
fprintf('θ(B,n): mean %.1f°, max %.1f°\n',mean(theta)*180/pi,max(theta)*180/pi);

%% -------------------- Helper --------------------
function A = ensure_RZ_order(Ain,nR,nZ)
    if isequal(size(Ain),[nR,nZ]), A=Ain;
    elseif isequal(size(Ain),[nZ,nR]), A=Ain';
    else, error('Array not [nR×nZ] or [nZ×nR]');
    end
end