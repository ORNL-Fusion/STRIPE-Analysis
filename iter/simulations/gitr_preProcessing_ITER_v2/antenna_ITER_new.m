%% ITER Antenna Surface Plasma Profiles + B-field (from profiles_iter_multi_GITRstyle_XY.nc)
% ======================================================================
% Robust against missing/invalid x,z axes (NetCDF fill-values).
% - Auto-detects coordinate variables if x/z are filled
% - Reorders all fields to match increasing R and Z
% - Uses griddedInterpolant on (R,Z)
% - Reduces patchiness with:
%     (1) optional light XY smoothing
%     (2) optional face-data smoothing on antenna mesh
%     (3) vertex-interpolated patch rendering
% ======================================================================

close all; clear; clc;

%% -------------------- User controls --------------------
neon_q          = 10;         % Neon charge state (0..10)
save_png        = true;      % save figures as 300-dpi PNGs
colormap_choice = turbo;     % colormap for plots
use_efit_B      = true;      % try EFIT B-field (guarded); else use nc Br/Bt/Bz

% ---- anti-patchiness controls ----
interp_method_xy   = 'linear';   % 'linear' safest; try 'makima' if supported by your MATLAB
doSmoothXY         = true;       % smooth XY fields before surface interpolation
sigma_xy           = 0.8;        % pixels on XY grid; try 0.6, 0.8, 1.0

doSmoothFaceData   = true;       % smooth antenna face data after centroid sampling
face_smooth_iters  = 2;          % 1-4 usually enough

use_vertex_interp  = true;       % render using vertex-interpolated colors
clip_negative_dens = true;       % clip tiny negative interpolation artifacts to zero

ncfile = 'profiles_iter_multi_GITRstyle_XY.nc';

%% -------------------- Load geometry --------------------
load('final_ITER_data.mat','planes','centroid','norm_vec');
fprintf('Loaded ITER geometry: %d triangular faces\n', size(planes,1));

r_centroid   = hypot(centroid(:,1), centroid(:,2));
phi_centroid = atan2(centroid(:,2), centroid(:,1));
z_centroid   = centroid(:,3);

%% -------------------- Load profiles --------------------
assert(isfile(ncfile), 'Missing NetCDF file: %s', ncfile);

% Read fields
ne    = ncread(ncfile,'ne');
te    = ncread(ncfile,'te');
ti    = ncread(ncfile,'ti');
br_nc = ncread(ncfile,'br');
bt_nc = ncread(ncfile,'bt');
bz_nc = ncread(ncfile,'bz');

niA = ncread(ncfile,'ni_all');
uRA = ncread(ncfile,'uR_all');
uZA = ncread(ncfile,'uZ_all');
uTA = ncread(ncfile,'uT_all');

Z_all = ncread(ncfile,'atomic_number');
q_all = ncread(ncfile,'charge_number');

%% -------------------- Determine grid sizes --------------------
sz_ne = size(ne);
if numel(sz_ne) ~= 2
    error('ne is not 2D. Got size=%s', mat2str(sz_ne));
end

ns = numel(Z_all);
s_ni = size(niA);
if numel(s_ni) ~= 3
    error('ni_all is not 3D. Got size=%s', mat2str(s_ni));
end

dims = s_ni;
if any(dims == ns)
    nzr = dims(dims ~= ns);
    if numel(nzr) ~= 2
        error('Cannot infer nR/nZ from ni_all dims=%s and ns=%d', mat2str(dims), ns);
    end
    nR_guess = nzr(1);
    nZ_guess = nzr(2);
else
    nR_guess = sz_ne(1);
    nZ_guess = sz_ne(2);
end

%% -------------------- Read coordinate axes --------------------
[R, Zg] = read_good_axes_from_nc(ncfile, nR_guess, nZ_guess);

fprintf('Loaded axes: nR=%d, nZ=%d\n', numel(R), numel(Zg));
fprintf('R(1)=%.6g  R(end)=%.6g | Z(1)=%.6g  Z(end)=%.6g\n', R(1), R(end), Zg(1), Zg(end));

nR = numel(R);
nZ = numel(Zg);

%% -------------------- Ensure ordering of fields --------------------
ne    = ensure_RZ_order(ne,    nR, nZ);
te    = ensure_RZ_order(te,    nR, nZ);
ti    = ensure_RZ_order(ti,    nR, nZ);
br_nc = ensure_RZ_order(br_nc, nR, nZ);
bt_nc = ensure_RZ_order(bt_nc, nR, nZ);
bz_nc = ensure_RZ_order(bz_nc, nR, nZ);

niA = ensure_RZS_order(niA, nR, nZ, ns);
uRA = ensure_RZS_order(uRA, nR, nZ, ns);
uZA = ensure_RZS_order(uZA, nR, nZ, ns);
uTA = ensure_RZS_order(uTA, nR, nZ, ns);

%% -------------------- Force increasing axes --------------------
[R, Zg, ne, te, ti, br_nc, bt_nc, bz_nc, niA, uRA, uZA, uTA] = ...
    force_increasing_axes(R, Zg, ne, te, ti, br_nc, bt_nc, bz_nc, niA, uRA, uZA, uTA);

%% -------------------- Optional XY smoothing --------------------
if doSmoothXY
    fprintf('Applying light XY smoothing (sigma = %.2f px)...\n', sigma_xy);

    ne    = local_imgaussfilt_fallback(ne,    sigma_xy);
    te    = local_imgaussfilt_fallback(te,    sigma_xy);
    ti    = local_imgaussfilt_fallback(ti,    sigma_xy);
    br_nc = local_imgaussfilt_fallback(br_nc, sigma_xy);
    bt_nc = local_imgaussfilt_fallback(bt_nc, sigma_xy);
    bz_nc = local_imgaussfilt_fallback(bz_nc, sigma_xy);

    for s = 1:ns
        niA(:,:,s) = local_imgaussfilt_fallback(niA(:,:,s), sigma_xy);
        uRA(:,:,s) = local_imgaussfilt_fallback(uRA(:,:,s), sigma_xy);
        uZA(:,:,s) = local_imgaussfilt_fallback(uZA(:,:,s), sigma_xy);
        uTA(:,:,s) = local_imgaussfilt_fallback(uTA(:,:,s), sigma_xy);
    end
end

%% -------------------- Quick sanity plot --------------------
fig0 = figure('Color','w','Name','Input ne and Antenna Centroids');
ax0  = axes(fig0);

imagesc(ax0, [R(1) R(end)], [Zg(1) Zg(end)], ne');
set(ax0,'YDir','normal');
hold(ax0,'on');
plot(ax0, r_centroid, z_centroid, 'k.', 'MarkerSize', 4);
hold(ax0,'off');
title(ax0,'Input n_e with Antenna Centroid Projection','FontSize',14,'FontWeight','bold');
xlabel(ax0,'R [m]'); ylabel(ax0,'Z [m]');
colorbar(ax0); colormap(ax0,colormap_choice);
axis(ax0,'tight');

%% -------------------- Species selection --------------------
idx_Dp  = find(Z_all==1  & q_all==1,      1, 'first');
idx_Neq = find(Z_all==10 & q_all==neon_q, 1, 'first');

bad_charge_meta = all(q_all(:)==Z_all(:));
if bad_charge_meta
    fprintf('Warning: charge_number metadata looks invalid (q==Z). Using index fallback ordering.\n');
end

if bad_charge_meta || isempty(idx_Dp) || isempty(idx_Neq)
    if ns >= 16
        idx_Dp  = 2;
        idx_Neq = 6 + neon_q;
    else
        error('Cannot infer species indices: ns=%d', ns);
    end
end
fprintf('Using species indices: D+ = %d, Ne^(q=%d) = %d\n', idx_Dp, neon_q, idx_Neq);

%% -------------------- Build interpolants --------------------
F_ne = griddedInterpolant({R, Zg}, ne, interp_method_xy, 'none');
F_te = griddedInterpolant({R, Zg}, te, interp_method_xy, 'none');
F_ti = griddedInterpolant({R, Zg}, ti, interp_method_xy, 'none');

F_br_nc = griddedInterpolant({R, Zg}, br_nc, interp_method_xy, 'none');
F_bt_nc = griddedInterpolant({R, Zg}, bt_nc, interp_method_xy, 'none');
F_bz_nc = griddedInterpolant({R, Zg}, bz_nc, interp_method_xy, 'none');

F_nD   = griddedInterpolant({R, Zg}, squeeze(niA(:,:,idx_Dp)),  interp_method_xy, 'none');
F_uRD  = griddedInterpolant({R, Zg}, squeeze(uRA(:,:,idx_Dp)),  interp_method_xy, 'none');
F_uZD  = griddedInterpolant({R, Zg}, squeeze(uZA(:,:,idx_Dp)),  interp_method_xy, 'none');
F_uTD  = griddedInterpolant({R, Zg}, squeeze(uTA(:,:,idx_Dp)),  interp_method_xy, 'none');

F_nNe  = griddedInterpolant({R, Zg}, squeeze(niA(:,:,idx_Neq)), interp_method_xy, 'none');
F_uRNe = griddedInterpolant({R, Zg}, squeeze(uRA(:,:,idx_Neq)), interp_method_xy, 'none');
F_uZNe = griddedInterpolant({R, Zg}, squeeze(uZA(:,:,idx_Neq)), interp_method_xy, 'none');
F_uTNe = griddedInterpolant({R, Zg}, squeeze(uTA(:,:,idx_Neq)), interp_method_xy, 'none');

%% -------------------- Interpolate to antenna centroids --------------------
ne_surf = F_ne(r_centroid, z_centroid);
te_surf = F_te(r_centroid, z_centroid);
ti_surf = F_ti(r_centroid, z_centroid);

%% -------------------- B-field interpolation (EFIT optional, guarded) --------------------
Br_surf = zeros(size(ne_surf));
Bt_surf = zeros(size(ne_surf));
Bz_surf = zeros(size(ne_surf));

if use_efit_B
    fprintf('Trying EFIT B-field from read_efit_data...\n');
    try
        read_efit_data;

        ok_vars = exist('Br','var') && exist('Bt','var') && exist('Bz','var') && ...
                  exist('r_efit','var') && exist('z_efit','var');

        if ok_vars
            FBr = griddedInterpolant({double(r_efit(:)), double(z_efit(:))}, double(Br').', 'linear', 'nearest');
            FBt = griddedInterpolant({double(r_efit(:)), double(z_efit(:))}, double(Bt').', 'linear', 'nearest');
            FBz = griddedInterpolant({double(r_efit(:)), double(z_efit(:))}, double(Bz').', 'linear', 'nearest');

            Br_surf = FBr(r_centroid, z_centroid);
            Bt_surf = FBt(r_centroid, z_centroid);
            Bz_surf = FBz(r_centroid, z_centroid);
            fprintf('EFIT B-field loaded.\n');
        else
            warning('EFIT Br/Bt/Bz variables not found; using NetCDF B-fields.');
            Br_surf = F_br_nc(r_centroid, z_centroid);
            Bt_surf = F_bt_nc(r_centroid, z_centroid);
            Bz_surf = F_bz_nc(r_centroid, z_centroid);
        end
    catch ME
        warning('EFIT B-field failed (%s). Using NetCDF B-fields.', ME.message);
        Br_surf = F_br_nc(r_centroid, z_centroid);
        Bt_surf = F_bt_nc(r_centroid, z_centroid);
        Bz_surf = F_bz_nc(r_centroid, z_centroid);
    end
else
    Br_surf = F_br_nc(r_centroid, z_centroid);
    Bt_surf = F_bt_nc(r_centroid, z_centroid);
    Bz_surf = F_bz_nc(r_centroid, z_centroid);
end
Bmag_surf = sqrt(Br_surf.^2 + Bt_surf.^2 + Bz_surf.^2);

%% -------------------- Species densities & velocities on surface --------------------
nD_surf   = F_nD(r_centroid, z_centroid);
uRD_surf  = F_uRD(r_centroid, z_centroid);
uZD_surf  = F_uZD(r_centroid, z_centroid);
uTD_surf  = F_uTD(r_centroid, z_centroid);
vDmag_surf = sqrt(uRD_surf.^2 + uZD_surf.^2 + uTD_surf.^2);

nNe_surf  = F_nNe(r_centroid, z_centroid);
uRNe_surf = F_uRNe(r_centroid, z_centroid);
uZNe_surf = F_uZNe(r_centroid, z_centroid);
uTNe_surf = F_uTNe(r_centroid, z_centroid);
vNeMag_surf = sqrt(uRNe_surf.^2 + uZNe_surf.^2 + uTNe_surf.^2);

%% -------------------- Sanitize --------------------
all_face_fields = {'ne_surf','te_surf','ti_surf','nD_surf','nNe_surf', ...
                   'uRD_surf','uZD_surf','uTD_surf','uRNe_surf','uZNe_surf','uTNe_surf', ...
                   'Br_surf','Bt_surf','Bz_surf','Bmag_surf','vDmag_surf','vNeMag_surf'};

for ii = 1:numel(all_face_fields)
    eval(sprintf('%s(~isfinite(%s)) = 0;', all_face_fields{ii}, all_face_fields{ii}));
end

if clip_negative_dens
    ne_surf = max(ne_surf, 0);
    te_surf = max(te_surf, 0);
    ti_surf = max(ti_surf, 0);
    nD_surf = max(nD_surf, 0);
    nNe_surf = max(nNe_surf, 0);
end

%% -------------------- Build patch mesh --------------------
nF = size(planes,1);
VA = planes(:,1:3);
VB = planes(:,4:6);
VC = planes(:,7:9);
V  = [VA; VB; VC];
F  = [(1:nF)' (1:nF)'+nF (1:nF)'+2*nF];

%% -------------------- Optional face smoothing on antenna mesh --------------------
if doSmoothFaceData
    fprintf('Applying face-data smoothing on antenna mesh (%d iterations)...\n', face_smooth_iters);

    ne_surf   = smooth_face_field(F, ne_surf,   face_smooth_iters);
    nD_surf   = smooth_face_field(F, nD_surf,   face_smooth_iters);
    nNe_surf  = smooth_face_field(F, nNe_surf,  face_smooth_iters);
    te_surf   = smooth_face_field(F, te_surf,   face_smooth_iters);
    ti_surf   = smooth_face_field(F, ti_surf,   face_smooth_iters);
    Bmag_surf = smooth_face_field(F, Bmag_surf, face_smooth_iters);

    Br_surf   = smooth_face_field(F, Br_surf,   face_smooth_iters);
    Bt_surf   = smooth_face_field(F, Bt_surf,   face_smooth_iters);
    Bz_surf   = smooth_face_field(F, Bz_surf,   face_smooth_iters);

    uRD_surf  = smooth_face_field(F, uRD_surf,  face_smooth_iters);
    uZD_surf  = smooth_face_field(F, uZD_surf,  face_smooth_iters);
    uTD_surf  = smooth_face_field(F, uTD_surf,  face_smooth_iters);

    uRNe_surf = smooth_face_field(F, uRNe_surf, face_smooth_iters);
    uZNe_surf = smooth_face_field(F, uZNe_surf, face_smooth_iters);
    uTNe_surf = smooth_face_field(F, uTNe_surf, face_smooth_iters);

    vDmag_surf  = sqrt(uRD_surf.^2 + uZD_surf.^2 + uTD_surf.^2);
    vNeMag_surf = sqrt(uRNe_surf.^2 + uZNe_surf.^2 + uTNe_surf.^2);
end

%% -------------------- Patch plots --------------------
Cset = {ne_surf, nD_surf, nNe_surf, te_surf, ti_surf, Bmag_surf};
titles = {'Electron Density n_e [m^{-3}]', ...
          'Deuterium Ion Density n_{D^+} [m^{-3}]', ...
          sprintf('Neon Density n_{Ne^{%d+}} [m^{-3}]', neon_q), ...
          'Electron Temperature T_e [eV]', ...
          'Ion Temperature T_i [eV]', ...
          'Magnetic Field Magnitude |B| [T]'};

for k = 1:numel(Cset)
    figk = figure('Color','w','Name',titles{k});
    axk  = axes(figk);

    if use_vertex_interp
        Cvertex = face_to_vertex_data(F, Cset{k}, size(V,1));
        patch(axk,'Faces',F,'Vertices',V,...
              'FaceVertexCData',Cvertex,...
              'FaceColor','interp','EdgeColor','none');
    else
        patch(axk,'Faces',F,'Vertices',V,...
              'FaceVertexCData',Cset{k},...
              'FaceColor','flat','EdgeColor','none');
    end

    axis(axk,'equal'); axis(axk,'tight'); axis(axk,'vis3d');
    xlabel(axk,'X [m]'); ylabel(axk,'Y [m]'); zlabel(axk,'Z [m]');
    title(axk,titles{k},'FontSize',15,'FontWeight','bold');
    colorbar(axk); colormap(axk,colormap_choice);
    view(axk,35,25);
    drawnow;

    if save_png
        exportgraphics(figk, sprintf('antenna_%02d.png', k), 'Resolution', 300);
    end
end

%% -------------------- B + velocity transforms, theta(B,n), fluxes --------------------
bx_surf = Br_surf .* cos(phi_centroid) - Bt_surf .* sin(phi_centroid);
by_surf = Br_surf .* sin(phi_centroid) + Bt_surf .* cos(phi_centroid);
bz_cart = Bz_surf;

b_mag = sqrt(bx_surf.^2 + by_surf.^2 + bz_cart.^2);
ubx = bx_surf ./ max(b_mag, 1e-30);
uby = by_surf ./ max(b_mag, 1e-30);
ubz = bz_cart ./ max(b_mag, 1e-30);

vx = uRD_surf .* cos(phi_centroid) - uTD_surf .* sin(phi_centroid);
vy = uRD_surf .* sin(phi_centroid) + uTD_surf .* cos(phi_centroid);
vz = uZD_surf;

vDmag = sqrt(vx.^2 + vy.^2 + vz.^2);
fluxD  = nD_surf  .* vDmag;
fluxNe = nNe_surf .* vNeMag_surf;

nmag = sqrt(sum(norm_vec.^2,2));
un = norm_vec ./ max(nmag, 1e-30);

dotBn = un(:,1).*ubx + un(:,2).*uby + un(:,3).*ubz;
dotBn = max(-1,min(1,dotBn));
theta = acos(dotBn);
theta(~isfinite(theta)) = 0;
theta(theta > pi/2) = pi - theta(theta > pi/2);

fprintf('Mean θ(B,n): %.2f°, Max θ(B,n): %.2f°\n', mean(theta)*180/pi, max(theta)*180/pi);

%% -------------------- Flux plots --------------------
flux_set = {fluxD, fluxNe};
flux_titles = { ...
    'D^+ Flux \Gamma_{D^+} = n_{D^+}|v_{D^+}| [m^{-2}s^{-1}]', ...
    sprintf('Ne^{%d+} Flux \\Gamma_{Ne} = n_{Ne}|v_{Ne}| [m^{-2}s^{-1}]', neon_q)};

for kf = 1:numel(flux_set)
    figf = figure('Color','w','Name',flux_titles{kf});
    axf  = axes(figf);

    if use_vertex_interp
        Cvertex = face_to_vertex_data(F, flux_set{kf}, size(V,1));
        patch(axf,'Faces',F,'Vertices',V,...
              'FaceVertexCData',Cvertex,...
              'FaceColor','interp','EdgeColor','none');
    else
        patch(axf,'Faces',F,'Vertices',V,...
              'FaceVertexCData',flux_set{kf},...
              'FaceColor','flat','EdgeColor','none');
    end

    axis(axf,'equal'); axis(axf,'tight'); axis(axf,'vis3d');
    xlabel(axf,'X [m]'); ylabel(axf,'Y [m]'); zlabel(axf,'Z [m]');
    title(axf,flux_titles{kf},'FontSize',15,'FontWeight','bold');
    colorbar(axf); colormap(axf,colormap_choice);
    view(axf,35,25);
    drawnow;

    if save_png
        exportgraphics(figf, sprintf('antenna_flux_%02d.png', kf), 'Resolution', 300);
    end
end

%% -------------------- Theta histogram --------------------
figh = figure('Color','w','Name','Theta(B,n) Histogram');
axh  = axes(figh);
histogram(axh, theta*180/pi);
title(axh,'\theta(B,n) Distribution on Antenna Surface','FontSize',14,'FontWeight','bold');
xlabel(axh,'\theta(B,n) [deg]'); ylabel(axh,'Count');

%% -------------------- Helper functions --------------------
function A = ensure_RZ_order(Ain,nR,nZ)
    if isequal(size(Ain),[nR,nZ])
        A = Ain;
    elseif isequal(size(Ain),[nZ,nR])
        A = Ain';
    else
        error('Array not [nR x nZ] or [nZ x nR]. Got %s', mat2str(size(Ain)));
    end
end

function A = ensure_RZS_order(Ain,nR,nZ,ns)
    s = size(Ain);
    if isequal(s,[nR,nZ,ns])
        A = Ain;
        return;
    end
    if numel(s) ~= 3
        error('Array not 3D for species field. Got %s', mat2str(s));
    end
    p = perms(1:3);
    for i = 1:size(p,1)
        if isequal(s(p(i,:)), [nR,nZ,ns])
            A = permute(Ain, p(i,:));
            return;
        end
    end
    error('Array not compatible with [nR x nZ x ns]. Got %s', mat2str(s));
end

function [R, Zg] = read_good_axes_from_nc(ncfile, nR_guess, nZ_guess)
    info = ncinfo(ncfile);
    vnames = string({info.Variables.Name});

    function ok = is_good_axis(vec, nExpected)
        vec = double(vec(:));
        ok = numel(vec)==nExpected && any(isfinite(vec)) && ...
             nnz(isfinite(vec))>2 && ...
             (max(vec(isfinite(vec))) < 1e30);
    end

    candR = ["x","r","R"];
    candZ = ["z","Z"];

    R = [];
    Zg = [];

    for nm = candR
        if any(vnames==nm)
            tmp = ncread(ncfile, char(nm));
            tmp = double(tmp(:));
            tmp(~isfinite(tmp) | tmp>1e30) = NaN;
            if is_good_axis(tmp, nR_guess)
                R = tmp;
                break;
            end
        end
    end

    for nm = candZ
        if any(vnames==nm)
            tmp = ncread(ncfile, char(nm));
            tmp = double(tmp(:));
            tmp(~isfinite(tmp) | tmp>1e30) = NaN;
            if is_good_axis(tmp, nZ_guess)
                Zg = tmp;
                break;
            end
        end
    end

    if isempty(R)
        for i = 1:numel(info.Variables)
            vv = info.Variables(i);
            if numel(vv.Size)==2 && any(vv.Size==1)
                n = max(vv.Size);
                if n == nR_guess
                    tmp = ncread(ncfile, vv.Name);
                    tmp = double(tmp(:));
                    tmp(~isfinite(tmp) | tmp>1e30) = NaN;
                    if is_good_axis(tmp, nR_guess)
                        fprintf('Axis R: using variable "%s"\n', vv.Name);
                        R = tmp;
                        break;
                    end
                end
            elseif numel(vv.Size)==1 && vv.Size==nR_guess
                tmp = ncread(ncfile, vv.Name);
                tmp = double(tmp(:));
                tmp(~isfinite(tmp) | tmp>1e30) = NaN;
                if is_good_axis(tmp, nR_guess)
                    fprintf('Axis R: using variable "%s"\n', vv.Name);
                    R = tmp;
                    break;
                end
            end
        end
    end

    if isempty(Zg)
        for i = 1:numel(info.Variables)
            vv = info.Variables(i);
            if numel(vv.Size)==2 && any(vv.Size==1)
                n = max(vv.Size);
                if n == nZ_guess
                    tmp = ncread(ncfile, vv.Name);
                    tmp = double(tmp(:));
                    tmp(~isfinite(tmp) | tmp>1e30) = NaN;
                    if is_good_axis(tmp, nZ_guess)
                        fprintf('Axis Z: using variable "%s"\n', vv.Name);
                        Zg = tmp;
                        break;
                    end
                end
            elseif numel(vv.Size)==1 && vv.Size==nZ_guess
                tmp = ncread(ncfile, vv.Name);
                tmp = double(tmp(:));
                tmp(~isfinite(tmp) | tmp>1e30) = NaN;
                if is_good_axis(tmp, nZ_guess)
                    fprintf('Axis Z: using variable "%s"\n', vv.Name);
                    Zg = tmp;
                    break;
                end
            end
        end
    end

    if isempty(R) || isempty(Zg)
        error(['Could not find valid coordinate vectors inside NetCDF. ' ...
               'Your export likely did not write x/z (they are all fill-values).']);
    end
end

function [R,Zg,ne,te,ti,br,bt,bz,niA,uRA,uZA,uTA] = ...
    force_increasing_axes(R,Zg,ne,te,ti,br,bt,bz,niA,uRA,uZA,uTA)

    R  = R(:);
    Zg = Zg(:);

    if R(end) <= R(1)
        fprintf('R axis descending; flipping R and all R-indexed fields.\n');
        R  = flipud(R);
        ne = flipud(ne); te = flipud(te); ti = flipud(ti);
        br = flipud(br); bt = flipud(bt); bz = flipud(bz);
        niA = flipud(niA); uRA = flipud(uRA); uZA = flipud(uZA); uTA = flipud(uTA);
    end

    if Zg(end) <= Zg(1)
        fprintf('Z axis descending; flipping Z and all Z-indexed fields.\n');
        Zg = flipud(Zg);
        ne = fliplr(ne); te = fliplr(te); ti = fliplr(ti);
        br = fliplr(br); bt = fliplr(bt); bz = fliplr(bz);
        niA = fliplr(niA); uRA = fliplr(uRA); uZA = fliplr(uZA); uTA = fliplr(uTA);
    end

    assert(R(end) > R(1),  'R not increasing after flip.');
    assert(Zg(end) > Zg(1), 'Z not increasing after flip.');
end

function A = local_imgaussfilt_fallback(Ain, sigma)
    if sigma <= 0
        A = Ain;
        return;
    end

    if exist('imgaussfilt','file') == 2
        A = imgaussfilt(Ain, sigma);
        return;
    end

    rad = max(1, ceil(3*sigma));
    x = -rad:rad;
    g = exp(-(x.^2)/(2*sigma^2));
    g = g / sum(g);
    K = g' * g;
    A = conv2(double(Ain), K, 'same');
end

function vtxData = face_to_vertex_data(F, faceData, nVertices)
    acc = zeros(nVertices,1);
    cnt = zeros(nVertices,1);

    for i = 1:size(F,1)
        vids = F(i,:);
        acc(vids) = acc(vids) + faceData(i);
        cnt(vids) = cnt(vids) + 1;
    end

    vtxData = acc ./ max(cnt,1);
end

function out = smooth_face_field(F, faceData, nIter)
    out = faceData(:);
    nF  = numel(out);

    adj = build_face_adjacency(F, nF);

    for it = 1:nIter
        newOut = out;
        for i = 1:nF
            nei = adj{i};
            if isempty(nei)
                continue;
            end
            newOut(i) = 0.5*out(i) + 0.5*mean(out(nei), 'omitnan');
        end
        out = newOut;
    end
end

function adj = build_face_adjacency(F, nF)
    adj = cell(nF,1);
    faceMap = containers.Map('KeyType','char','ValueType','any');

    for i = 1:nF
        tri = F(i,:);
        edges = [tri([1 2]); tri([2 3]); tri([3 1])];
        for e = 1:3
            a = min(edges(e,:));
            b = max(edges(e,:));
            key = sprintf('%d_%d', a, b);
            if isKey(faceMap, key)
                faceMap(key) = [faceMap(key), i];
            else
                faceMap(key) = i;
            end
        end
    end

    keysList = keys(faceMap);
    for k = 1:numel(keysList)
        facesHere = faceMap(keysList{k});
        facesHere = unique(facesHere);
        if numel(facesHere) >= 2
            for i = 1:numel(facesHere)
                fi = facesHere(i);
                others = facesHere(facesHere ~= fi);
                adj{fi} = unique([adj{fi}, others]);
            end
        end
    end
end