%% ITER Antenna Surface Plasma Profiles + B-field (from profiles.nc)
% ======================================================================
% Uses a single NetCDF file (profiles.nc) containing:
%   x,z,ne,te,ti,br,bt,bz,
%   ni_all(R,Z,S), uR_all(R,Z,S), uZ_all(R,Z,S), uT_all(R,Z,S),
%   atomic_number(S), charge_number(S)
% Maps data to ITER antenna triangles and makes high-resolution patch plots.
% ======================================================================

% clear; clc; close all;

%% -------------------- User controls --------------------
neon_q  = 9;            % Neon charge state (0..10)
save_png = true;         % save figures as 300-dpi PNGs
colormap_choice = turbo; % colormap used for all plots
use_efit_B = true;       % true: use read_efit_data for Br/Bt/Bz on surface

%% -------------------- Load geometry --------------------
load('final_ITER_data.mat','planes','centroid','norm_vec','X','Y','Z');
fprintf('Loaded ITER geometry: %d triangular faces\n', size(planes,1));

r_centroid = hypot(centroid(:,1),centroid(:,2));
phi_centroid = atan2(centroid(:,2),centroid(:,1));

%% -------------------- Load profiles --------------------
ncfile = 'profiles_iter_multi_GITRstyle_XY.nc';

% assert(isfile(ncfile),'profiles_iter_wideGrid.nc not found');

R  = ncread(ncfile,'x');           % [nR]
Zg = ncread(ncfile,'z');           % [nZ]
[Rg, Zm] = meshgrid(R, Zg);

ne   = ncread(ncfile,'ne');  te = ncread(ncfile,'te');  ti = ncread(ncfile,'ti');
br_nc = ncread(ncfile,'br'); bt_nc = ncread(ncfile,'bt'); bz_nc = ncread(ncfile,'bz');

niA  = ncread(ncfile,'ni_all');
uRA  = ncread(ncfile,'uR_all');
uZA  = ncread(ncfile,'uZ_all');
uTA  = ncread(ncfile,'uT_all');

Z_all = ncread(ncfile,'atomic_number');
q_all = ncread(ncfile,'charge_number');
[nR,nZ,ns] = size(niA);
fprintf('SOLPS grid: nR=%d, nZ=%d, ns=%d\n', nR,nZ,ns);

% Ensure all 2D fields are in [nR x nZ]
ne = ensure_RZ_order(ne, numel(R), numel(Zg));
te = ensure_RZ_order(te, numel(R), numel(Zg));
ti = ensure_RZ_order(ti, numel(R), numel(Zg));
br_nc = ensure_RZ_order(br_nc, numel(R), numel(Zg));
bt_nc = ensure_RZ_order(bt_nc, numel(R), numel(Zg));
bz_nc = ensure_RZ_order(bz_nc, numel(R), numel(Zg));
niA = ensure_RZS_order(niA, numel(R), numel(Zg), ns);
uRA = ensure_RZS_order(uRA, numel(R), numel(Zg), ns);
uZA = ensure_RZS_order(uZA, numel(R), numel(Zg), ns);
uTA = ensure_RZS_order(uTA, numel(R), numel(Zg), ns);

fig0 = figure('Color','w','Name','Input ne and Antenna Centroids');
ax0 = axes(fig0);
imagesc(ax0, R, Zg, ne'); hold(ax0,'on');
plot(r_centroid, centroid(:,3), 'k.');
set(ax0,'YDir','normal');
title(ax0,'Input n_e with Antenna Centroid Projection','FontSize',14,'FontWeight','bold');
xlabel('R [m]');
ylabel('Z [m]');
colorbar;

swap = @(A) permute(A,[2 1]); % helper for interp2(Z,R)

%% -------------------- Species selection --------------------
% Prefer metadata-based lookup, but fall back if charge metadata is invalid.
idx_Dp  = find(Z_all==1 & q_all==1, 1, 'first');      % D+
idx_Neq = find(Z_all==10 & q_all==neon_q, 1, 'first');% Ne^q+

% Detect common bad metadata case: charge_number copied from atomic_number.
bad_charge_meta = all(q_all(:)==Z_all(:));
if bad_charge_meta
    fprintf('Warning: charge_number metadata appears invalid (q==Z for all species). Using index fallback.\n');
end

if bad_charge_meta || isempty(idx_Dp) || isempty(idx_Neq)
    % Expected ordering in this workflow: [D0 D+ He0 He1 He2 Ne0..Ne10]
    if ns >= 16
        idx_Dp = 2;
        idx_Neq = 6 + neon_q;   % Ne0 at 6, Ne10 at 16
    else
        error('Cannot infer species indices: ns=%d is not compatible with fallback ordering.', ns);
    end
end

if idx_Dp < 1 || idx_Dp > ns || idx_Neq < 1 || idx_Neq > ns
    error('Selected species indices out of bounds: idx_Dp=%d idx_Neq=%d ns=%d', idx_Dp, idx_Neq, ns);
end
fprintf('Using species: D+ = %d, Ne^(q=%d) = %d\n',idx_Dp,neon_q,idx_Neq);
fprintf('Selected species metadata: D+(Z=%g,q=%g), Ne(Z=%g,q=%g)\n', ...
        Z_all(idx_Dp), q_all(idx_Dp), Z_all(idx_Neq), q_all(idx_Neq));

%% -------------------- Interpolate to antenna centroids --------------------
method = 'natural';

ne_surf = interp2(Rg,Zm,swap(ne), r_centroid,centroid(:,3),method);
te_surf = interp2(Rg,Zm,swap(te), r_centroid,centroid(:,3),method);
ti_surf = interp2(Rg,Zm,swap(ti), r_centroid,centroid(:,3),method);

if use_efit_B
    fprintf('Using EFIT B-field from read_efit_data for surface interpolation...\n');
    read_efit_data; % provides Br, Bt, Bz, r_efit, z_efit
    FBr = griddedInterpolant({z_efit, r_efit}, Br, 'linear', 'nearest');
    FBt = griddedInterpolant({z_efit, r_efit}, Bt, 'linear', 'nearest');
    FBz = griddedInterpolant({z_efit, r_efit}, Bz, 'linear', 'nearest');
    Br_surf = FBr(centroid(:,3), r_centroid);
    Bt_surf = FBt(centroid(:,3), r_centroid);
    Bz_surf = FBz(centroid(:,3), r_centroid);
else
    fprintf('Using Br/Bt/Bz from %s\n', ncfile);
    Br_surf = interp2(Rg,Zm,swap(br_nc), r_centroid,centroid(:,3),method);
    Bt_surf = interp2(Rg,Zm,swap(bt_nc), r_centroid,centroid(:,3),method);
    Bz_surf = interp2(Rg,Zm,swap(bz_nc), r_centroid,centroid(:,3),method);
end
Bmag_surf = sqrt(Br_surf.^2 + Bt_surf.^2 + Bz_surf.^2);

nD_surf  = interp2(Rg,Zm,swap(niA(:,:,idx_Dp)),  r_centroid,centroid(:,3),method);
uRD_surf = interp2(Rg,Zm,swap(uRA(:,:,idx_Dp)),  r_centroid,centroid(:,3),method);
uZD_surf = interp2(Rg,Zm,swap(uZA(:,:,idx_Dp)),  r_centroid,centroid(:,3),method);
uTD_surf = interp2(Rg,Zm,swap(uTA(:,:,idx_Dp)),  r_centroid,centroid(:,3),method);
vDmag_surf = sqrt(uRD_surf.^2 + uZD_surf.^2 + uTD_surf.^2);

nNe_surf  = interp2(Rg,Zm,swap(niA(:,:,idx_Neq)), r_centroid,centroid(:,3),method);
uRNe_surf = interp2(Rg,Zm,swap(uRA(:,:,idx_Neq)), r_centroid,centroid(:,3),method);
uZNe_surf = interp2(Rg,Zm,swap(uZA(:,:,idx_Neq)), r_centroid,centroid(:,3),method);
uTNe_surf = interp2(Rg,Zm,swap(uTA(:,:,idx_Neq)), r_centroid,centroid(:,3),method);
vNeMag_surf = sqrt(uRNe_surf.^2 + uZNe_surf.^2 + uTNe_surf.^2);

% sanitize non-finite interpolated values
ne_surf(~isfinite(ne_surf)) = 0; te_surf(~isfinite(te_surf)) = 0; ti_surf(~isfinite(ti_surf)) = 0;
nD_surf(~isfinite(nD_surf)) = 0; nNe_surf(~isfinite(nNe_surf)) = 0;
uRD_surf(~isfinite(uRD_surf)) = 0; uZD_surf(~isfinite(uZD_surf)) = 0; uTD_surf(~isfinite(uTD_surf)) = 0;
uRNe_surf(~isfinite(uRNe_surf)) = 0; uZNe_surf(~isfinite(uZNe_surf)) = 0; uTNe_surf(~isfinite(uTNe_surf)) = 0;
Br_surf(~isfinite(Br_surf)) = 0; Bt_surf(~isfinite(Bt_surf)) = 0; Bz_surf(~isfinite(Bz_surf)) = 0; Bmag_surf(~isfinite(Bmag_surf)) = 0;
vDmag_surf(~isfinite(vDmag_surf)) = 0; vNeMag_surf(~isfinite(vNeMag_surf)) = 0;

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
Cset = {ne_surf, nD_surf, nNe_surf, te_surf, ti_surf, Bmag_surf};
titles = {'Electron Density n_e [m^{-3}]', ...
          'Deuterium Ion Density n_{D^+} [m^{-3}]', ...
          sprintf('Neon Density n_{Ne^{%d+}} [m^{-3}]', neon_q), ...
          'Electron Temperature T_e [eV]', ...
          'Ion Temperature T_i [eV]', ...
          'Magnetic Field Magnitude |B| [T]'};

for k = 1:numel(Cset)
    figk = figure('Color','w','Name',titles{k});
    axk = axes(figk);
    patch(axk,'Faces',F,'Vertices',V,...
          'FaceVertexCData',Cset{k},...   % one color per face
          'FaceColor','flat','EdgeColor','none');
    axis(axk,'equal'); axis(axk,'tight'); axis(axk,'vis3d');
    xlabel(axk,'X [m]'); ylabel(axk,'Y [m]'); zlabel(axk,'Z [m]');
    title(axk,titles{k},'FontSize',15,'FontWeight','bold');
    colorbar(axk); colormap(axk,colormap_choice);
    view(axk,35,25);
    drawnow;
    if save_png
        outname = sprintf('antenna_%02d.png', k);
        exportgraphics(gcf, outname, 'Resolution', 300);
    end
end


% %% -------------------- Diagnostics --------------------
% fprintf('\n=== ITER Antenna Diagnostics ===\n');
% fprintf('n_e: %.2e–%.2e [m^-3]\n',min(ne_surf),max(ne_surf));
% fprintf('n_D^+: %.2e–%.2e [m^-3]\n',min(ni_D_surf),max(ni_D_surf));
% fprintf('n_Ne^0: %.2e–%.2e [m^-3]\n',min(ni_Ne0_surf),max(ni_Ne0_surf));
% fprintf('n_Ne^{10+}: %.2e–%.2e [m^-3]\n',min(ni_Ne10_surf),max(ni_Ne10_surf));
% fprintf('θ(B,n): mean %.1f°, max %.1f°\n',mean(theta)*180/pi,max(theta)*180/pi);
%% -------------------- Magnetic field & velocity transformations --------------------
fprintf('Performing magnetic and velocity transformations...\n');

% --- Magnetic field transformation (R,Z,phi → X,Y,Z) ---
phi_centroid = atan2(centroid(:,2), centroid(:,1));

bx_surf = Br_surf .* cos(phi_centroid) - Bt_surf .* sin(phi_centroid);
by_surf = Br_surf .* sin(phi_centroid) + Bt_surf .* cos(phi_centroid);
bz_surf = Bz_surf;

b_mag_surf = sqrt(bx_surf.^2 + by_surf.^2 + bz_surf.^2);
ubx_surf = bx_surf ./ max(b_mag_surf, 1e-30);
uby_surf = by_surf ./ max(b_mag_surf, 1e-30);
ubz_surf = bz_surf ./ max(b_mag_surf, 1e-30);

% --- Velocity transformation (R,Z,phi → X,Y,Z) ---
% If you already interpolated vr_surf, vt_surf, vz_surf
% (add them to your surface interpolation block)
vx_surf = uRD_surf .* cos(phi_centroid) - uTD_surf .* sin(phi_centroid);
vy_surf = uRD_surf .* sin(phi_centroid) + uTD_surf .* cos(phi_centroid);
vz_surf = uZD_surf;  % keep Z-component as-is

v_mag_surf = sqrt(vx_surf.^2 + vy_surf.^2 + vz_surf.^2);
uvx_surf = vx_surf ./ max(v_mag_surf, 1e-30);
uvy_surf = vy_surf ./ max(v_mag_surf, 1e-30);
uvz_surf = vz_surf ./ max(v_mag_surf, 1e-30);

% --- Example flux (O8+ shown generically, adjust species name as needed) ---
fluxD_surf  = nD_surf  .* vDmag_surf;   % D+ particle flux [m^-2 s^-1]
fluxNe_surf = nNe_surf .* vNeMag_surf;  % Ne^q+ particle flux [m^-2 s^-1]

% --- Surface normals and θ(B,n) calculation ---
norm_vec_mag = sqrt(norm_vec(:,1).^2 + norm_vec(:,2).^2 + norm_vec(:,3).^2);
unorm_vec = norm_vec ./ norm_vec_mag;

theta = acos(unorm_vec(:,1).*ubx_surf + ...
                unorm_vec(:,2).*uby_surf + ...
                unorm_vec(:,3).*ubz_surf);
theta(isnan(theta)) = 0;
theta(theta > pi/2) = pi - theta(theta > pi/2);

fprintf('Mean θ(B,n): %.2f°, Max θ(B,n): %.2f°\n', ...
        mean(theta)*180/pi, max(theta)*180/pi);

%% -------------------- Flux plots on antenna surface --------------------
flux_set = {fluxD_surf, fluxNe_surf};
flux_titles = { ...
    'D^+ Flux \Gamma_{D^+} = n_{D^+}|v_{D^+}| [m^{-2}s^{-1}]', ...
    sprintf('Ne^{%d+} Flux \\Gamma_{Ne} = n_{Ne}|v_{Ne}| [m^{-2}s^{-1}]', neon_q)};

for kf = 1:numel(flux_set)
    figf = figure('Color','w','Name',flux_titles{kf});
    axf = axes(figf);
    patch(axf,'Faces',F,'Vertices',V,...
          'FaceVertexCData',flux_set{kf},...
          'FaceColor','flat','EdgeColor','none');
    axis(axf,'equal'); axis(axf,'tight'); axis(axf,'vis3d');
    xlabel(axf,'X [m]'); ylabel(axf,'Y [m]'); zlabel(axf,'Z [m]');
    title(axf,flux_titles{kf},'FontSize',15,'FontWeight','bold');
    colorbar(axf); colormap(axf,colormap_choice);
    view(axf,35,25);
    drawnow;
    if save_png
        outname = sprintf('antenna_flux_%02d.png', kf);
        exportgraphics(figf, outname, 'Resolution', 300);
    end
end

%% --- Visualization ---
% Normal and B-field vectors
figv = figure('Color','w','Name','Surface Normals and B Vectors');
axv = axes(figv);
quiver3(axv,centroid(:,1), centroid(:,2), centroid(:,3), norm_vec(:,1), norm_vec(:,2), norm_vec(:,3)); 
hold(axv,'on');
quiver3(axv,centroid(:,1), centroid(:,2), centroid(:,3), bx_surf, by_surf, bz_surf);
title(axv,'Surface Normals and Magnetic Field Vectors','FontSize',14,'FontWeight','bold');
xlabel(axv,'X [m]'); ylabel(axv,'Y [m]'); zlabel(axv,'Z [m]');
axis(axv,'equal'); axis(axv,'tight'); axis(axv,'vis3d');
legend(axv,{'Surface normal','B vector'}, 'Location', 'best');
drawnow;

% Histogram of theta
figh = figure('Color','w','Name','Theta(B,n) Histogram');
axh = axes(figh);
histogram(axh,theta*180/pi);
title(axh,'\theta(B,n) Distribution on Antenna Surface','FontSize',14,'FontWeight','bold');
xlabel(axh,'\theta(B,n) [deg]');
ylabel(axh,'Count');
drawnow;
%% -------------------- Helper --------------------
function A = ensure_RZ_order(Ain,nR,nZ)
    if isequal(size(Ain),[nR,nZ]), A=Ain;
    elseif isequal(size(Ain),[nZ,nR]), A=Ain';
    else, error('Array not [nR×nZ] or [nZ×nR]');
    end
end

function A = ensure_RZS_order(Ain,nR,nZ,ns)
    s = size(Ain);
    if isequal(s,[nR,nZ,ns])
        A = Ain;
        return;
    end
    if numel(s)~=3
        error('Array not 3D for species field');
    end
    p = perms(1:3);
    for i = 1:size(p,1)
        if isequal(s(p(i,:)), [nR,nZ,ns])
            A = permute(Ain, p(i,:));
            return;
        end
    end
    error('Array not compatible with [nR x nZ x ns]');
end
