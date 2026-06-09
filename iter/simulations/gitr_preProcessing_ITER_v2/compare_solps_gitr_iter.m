%% compare_SOLPSmat_vs_GITRstyle_midplane.m
% Compare midplane profiles from SOLPS .mat (painted to rect grid)
% vs an external GITR-style NetCDF (x,z,ne,ni_all).
% Outputs 3 overlay plots: ne, nD+, nNe_total at Z~0.

clear; clc; close all;

%% ------------------- Inputs -------------------
matFile     = 'solps_iter.mat';
gitr_cmp_nc = 'profiles_iter_multi_GITRstyle_XY.nc';

% Rectangular grid for painting SOLPS cells
nR = 800;
nZ = 400;

% Species indices (your convention)
iDplus = 2;      % D+
idxNe  = 6:16;   % Ne0..Ne10+

%% ------------------- Load SOLPS .mat -------------------
S = load(matFile);
Geo   = S.Geo;
State = S.State;

r = double(Geo.pr);  % [4 x Nc]
z = double(Geo.pz);  % [4 x Nc]
Nc = size(r,2);

ne_cv = double(State.ne(:));    % [Nc x 1]
na_cv = double(State.na);       % [Nc x ns]
ns    = size(na_cv,2);

fprintf('Loaded SOLPS: Nc=%d cells, ns=%d species from %s\n', Nc, ns, matFile);

%% ------------------- Define SOLPS painted (R,Z) grid -------------------
rmin = min(r(:)); rmax = max(r(:));
zmin = min(z(:)); zmax = max(z(:));

rgrid = linspace(rmin, rmax, nR);
zgrid = linspace(zmin, zmax, nZ);
[Rg, Zg] = meshgrid(rgrid, zgrid);   % [nZ x nR]

%% ------------------- Paint helper -------------------
fill_by_cells = @(vals) paintCells(vals, r, z, Rg, Zg);

%% ------------------- Paint SOLPS fields -------------------
ne_q = fill_by_cells(ne_cv);         % [nZ x nR]

ni_q = cell(1,ns);
for s = 1:ns
    ni_q{s} = fill_by_cells(na_cv(:,s));
end

% Clean
ne_q(ne_q<=0) = NaN;
for s = 1:ns
    tmp = ni_q{s}; tmp(tmp<=0) = NaN; ni_q{s} = tmp;
end

%% ------------------- Midplane slice from SOLPS painted grid -------------------
[~, iz_SOLPS] = min(abs(zgrid - 0));
fprintf('SOLPS midplane: iz=%d (Z=%.4f m)\n', iz_SOLPS, zgrid(iz_SOLPS));

ne_SOLPS_mid = ne_q(iz_SOLPS,:);

if iDplus <= ns
    nD_SOLPS_mid = ni_q{iDplus}(iz_SOLPS,:);
else
    nD_SOLPS_mid = NaN(size(ne_SOLPS_mid));
    warning('iDplus=%d exceeds ns=%d in SOLPS.', iDplus, ns);
end

nNe_SOLPS_mid = zeros(size(ne_SOLPS_mid));
for k = idxNe
    if k <= ns
        tmp = ni_q{k}(iz_SOLPS,:);
        tmp(~isfinite(tmp)) = 0;
        nNe_SOLPS_mid = nNe_SOLPS_mid + tmp;
    end
end
nNe_SOLPS_mid(nNe_SOLPS_mid<=0) = NaN;

%% ------------------- Read EXTERNAL GITR-style NetCDF -------------------
assert(exist(gitr_cmp_nc,'file')==2, 'External GITR file not found: %s', gitr_cmp_nc);

xG = ncread(gitr_cmp_nc,'x');   % [nX]
zG = ncread(gitr_cmp_nc,'z');   % [nZ]
nXg = numel(xG); nZg = numel(zG);

neG_raw = ncread(gitr_cmp_nc,'ne');     % typical [nX x nZ] or [nZ x nX]
neG = local_force_ZxX(neG_raw, nXg, nZg);   % [nZ x nX]

[~, iz_G] = min(abs(zG - 0));
fprintf('GITR midplane: iz=%d (Z=%.4f m)\n', iz_G, zG(iz_G));

ne_G_mid = neG(iz_G,:); ne_G_mid(ne_G_mid<=0) = NaN;

% ni_all if present
vars = string({ncinfo(gitr_cmp_nc).Variables.Name});
hasNi = any(vars=="ni_all");

nD_G_mid  = [];
nNe_G_mid = [];
if hasNi
    niG_all = ncread(gitr_cmp_nc,'ni_all'); % [nX x nZ x ns] or [nZ x nX x ns] etc

    if size(niG_all,3) >= iDplus
        nD2D = local_force_ZxX_slice(niG_all, nXg, nZg, iDplus);   % [nZ x nX]
        nD_G_mid = nD2D(iz_G,:); nD_G_mid(nD_G_mid<=0)=NaN;
    else
        warning('GITR ni_all has only %d species; iDplus=%d invalid.', size(niG_all,3), iDplus);
    end

    if size(niG_all,3) >= max(idxNe)
        nNe2D = zeros(nZg, nXg);
        for k = idxNe
            nNe2D = nNe2D + local_force_ZxX_slice(niG_all, nXg, nZg, k);
        end
        nNe_G_mid = nNe2D(iz_G,:); nNe_G_mid(nNe_G_mid<=0)=NaN;
    else
        warning('GITR ni_all has only %d species; cannot sum idxNe up to %d.', size(niG_all,3), max(idxNe));
    end
else
    warning('External GITR file has no ni_all; only n_e will be compared.');
end

%% ------------------- Plot overlays -------------------
figure('Color','w','Position',[120 120 1350 420]);
tiledlayout(1,3,'TileSpacing','compact','Padding','compact');

nexttile;
semilogy(rgrid, max(ne_SOLPS_mid,1), 'k-',  'LineWidth', 2); hold on; grid on; box on;
semilogy(xG,    max(ne_G_mid,1),     'k--', 'LineWidth', 2);
xlabel('R [m]'); ylabel('n_e [m^{-3}]'); title('Midplane n_e (SOLPS .mat vs GITR)');
legend('SOLPS(.mat painted)','GITR(NetCDF)','Location','best');

nexttile;
if ~isempty(nD_G_mid)
    semilogy(rgrid, max(nD_SOLPS_mid,1), 'b-',  'LineWidth', 2); hold on; grid on; box on;
    semilogy(xG,    max(nD_G_mid,1),     'b--', 'LineWidth', 2);
    xlabel('R [m]'); ylabel('n_{D^+} [m^{-3}]'); title('Midplane n_{D^+}');
    legend('SOLPS','GITR','Location','best');
else
    axis off; text(0.1,0.5,'No D+ available in GITR ni\_all', 'FontSize', 13);
end

nexttile;
if ~isempty(nNe_G_mid)
    semilogy(rgrid, max(nNe_SOLPS_mid,1), 'r-',  'LineWidth', 2); hold on; grid on; box on;
    semilogy(xG,    max(nNe_G_mid,1),     'r--', 'LineWidth', 2);
    xlabel('R [m]'); ylabel('n_{Ne,total} [m^{-3}]'); title('Midplane n_{Ne,total}');
    legend('SOLPS','GITR','Location','best');
else
    axis off; text(0.1,0.5,'No Ne-total available in GITR ni\_all', 'FontSize', 13);
end

sgtitle('Midplane compare at Z \approx 0: SOLPS (.mat painted) vs external GITR-style NetCDF','FontSize',14);

fprintf('✅ Done: SOLPS(.mat) vs GITR midplane overlays.\n');
%% =============================================================================
%% === ANTENNA MIDPLANE: Neon charge-state abundance histogram (using centroids) ===
doAntennaMidplaneNeHist = true;

geom_mat = 'final_ITER_data.mat';   % provides planes, centroid, norm_vec

% --- Selection knobs (edit as needed) ---
z_mid_tol = 0.02;          % [m] keep antenna triangles with |Z| <= z_mid_tol  (midplane band)
useRwindow = true;
R_ant_min  = 8.20;         % [m] optional antenna R window
R_ant_max  = 8.36;         % [m]

% If you want to also filter by phi, enable and set window:
usePhiWindow = false;
phi_min = -pi;             % [rad]
phi_max =  pi;             % [rad]

NeNames = {'Ne^{0}','Ne^{1+}','Ne^{2+}','Ne^{3+}','Ne^{4+}', ...
           'Ne^{5+}','Ne^{6+}','Ne^{7+}','Ne^{8+}','Ne^{9+}','Ne^{10+}'};

if doAntennaMidplaneNeHist && exist(geom_mat,'file')
    load(geom_mat,'planes','centroid','norm_vec'); %#ok<LOAD>
    fprintf('\nLoaded ITER geometry: %d triangular faces\n', size(planes,1));

    % Cylindrical coords of triangle centroids
    r_centroid   = hypot(centroid(:,1), centroid(:,2));
    phi_centroid = atan2(centroid(:,2), centroid(:,1));
    z_centroid   = centroid(:,3);

    % --- Antenna midplane selection mask ---
    m = abs(z_centroid) <= z_mid_tol;
    if useRwindow
        m = m & (r_centroid >= R_ant_min) & (r_centroid <= R_ant_max);
    end
    if usePhiWindow
        m = m & (phi_centroid >= phi_min) & (phi_centroid <= phi_max);
    end

    fprintf('Antenna-midplane centroid selection: kept %d / %d triangles (|Z|<=%.3f m)\n', ...
        nnz(m), numel(m), z_mid_tol);

    if nnz(m) < 5
        warning('Too few antenna-midplane triangles selected. Loosen z_mid_tol / windows.');
    else
        Rq = r_centroid(m);
        Zq = z_centroid(m);

        % --- Interpolants from SOLPS painted grid ---
        % ne_q, ni_q{species} are on (zgrid,rgrid) with size [nZ x nR]
        % IMPORTANT: griddedInterpolant expects axes as {zgrid,rgrid} for arrays [nZ x nR]
        Fne = griddedInterpolant({zgrid, rgrid}, ne_q, 'linear', 'none');
        ne_ant = Fne(Zq, Rq);

        % Neon densities at antenna midplane triangles
        nNe_q_ant = NaN(numel(Rq), numel(idxNe));
        for j = 1:numel(idxNe)
            k = idxNe(j);
            if k <= ns
                Fni = griddedInterpolant({zgrid, rgrid}, ni_q{k}, 'linear', 'none');
                nNe_q_ant(:,j) = Fni(Zq, Rq);
            else
                nNe_q_ant(:,j) = 0;
            end
        end

        % Total Ne and fractions per triangle
        nNe_tot = sum(nNe_q_ant, 2, 'omitnan');
        f_q = nNe_q_ant ./ max(nNe_tot, 1e-60);
        f_q(~isfinite(f_q)) = 0;

        % Triangle-averaged abundance (simple mean of f_q)
        fq_mean = mean(f_q, 1, 'omitnan');
        fq_mean = fq_mean ./ max(sum(fq_mean), 1e-60);

        fprintf('\n=== Antenna-midplane Neon charge-state abundance (triangle-averaged) ===\n');
        if useRwindow
            fprintf('R window: [%.3f, %.3f] m\n', R_ant_min, R_ant_max);
        end
        if usePhiWindow
            fprintf('phi window: [%.3f, %.3f] rad\n', phi_min, phi_max);
        end
        fprintf('midplane band: |Z| <= %.3f m\n', z_mid_tol);

        for j = 1:numel(idxNe)
            q = j-1; % Ne0..Ne10+
            fprintf('%4s (q=%2d): %8.4f  (%6.2f%%)\n', NeNames{j}, q, fq_mean(j), 100*fq_mean(j));
        end

        % --- HISTOGRAM-style plot (bar chart over q=0..10) ---
        figure('Color','w','Position',[140 140 900 450]);
        bar(0:10, 100*fq_mean, 'FaceColor','flat');
        xlabel('Ne charge state q (0..10)');
        ylabel('Abundance [%]');
        title(sprintf('Antenna-midplane Neon charge-state abundance (|Z|<=%.2f m)', z_mid_tol));
        grid on; box on;

        % Optional: show points used
        figure('Color','w','Position',[160 160 1050 420]);
        tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

        nexttile;
        scatter(Rq, Zq, 8, log10(max(nNe_tot,1)), 'filled'); colorbar;
        xlabel('R [m]'); ylabel('Z [m]'); axis equal tight;
        title('Selected antenna-midplane triangles: log_{10}(n_{Ne,total})');

        nexttile;
        scatter(Rq, Zq, 8, 100*(nNe_tot ./ max(ne_ant,1e-60)), 'filled'); colorbar;
        xlabel('R [m]'); ylabel('Z [m]'); axis equal tight;
        title('Selected triangles: 100×(n_{Ne,total}/n_e)');
    end
else
    if doAntennaMidplaneNeHist
        warning('Antenna-midplane Neon histogram skipped: %s not found.', geom_mat);
    end
end
%% =============================================================================
%% ------------------- Local helpers -------------------
function grid_vals = paintCells(cv_vals, r4xN, z4xN, X, Z)
% Paint SOLPS quad cells to a rect grid with constant per-cell values.
    grid_vals = NaN(size(X));
    N = size(r4xN,2);
    for i = 1:N
        ri = r4xN(:,i);
        zi = z4xN(:,i);
        if any(~isfinite(ri)) || any(~isfinite(zi)), continue; end
        [in,on] = inpolygon(X, Z, ri, zi);
        m = (in | on);
        if any(m(:))
            grid_vals(m) = cv_vals(i);
        end
    end
end

function V = local_force_ZxX(A, nX, nZ)
% Coerce A to [nZ x nX]
    sz = size(A);
    if isequal(sz, [nZ nX])
        V = A;
    elseif isequal(sz, [nX nZ])
        V = A.';
    elseif numel(A) == nX*nZ
        V = reshape(A, [nX nZ]).';
    else
        error('local_force_ZxX: cannot coerce size %s to [nZ x nX]=[%d %d].', mat2str(sz), nZ, nX);
    end
end

function Vk = local_force_ZxX_slice(A3, nX, nZ, k)
% Extract species slice and coerce to [nZ x nX]
    sz = size(A3);
    if numel(sz) ~= 3
        error('local_force_ZxX_slice: expected 3D array, got %s', mat2str(sz));
    end

    if sz(1)==nX && sz(2)==nZ
        Vk = A3(:,:,k).';
        return
    elseif sz(1)==nZ && sz(2)==nX
        Vk = A3(:,:,k);
        return
    end

    perms_all = perms(1:3);
    for i = 1:size(perms_all,1)
        P = perms_all(i,:);
        B = permute(A3, P);
        s2 = size(B);
        if s2(1)==nX && s2(2)==nZ
            Vk = B(:,:,k).';
            return
        elseif s2(1)==nZ && s2(2)==nX
            Vk = B(:,:,k);
            return
        end
    end

    error('local_force_ZxX_slice: cannot coerce size %s with nX=%d nZ=%d.', mat2str(sz), nX, nZ);
end