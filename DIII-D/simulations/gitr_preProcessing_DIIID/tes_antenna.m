close all;
clear all;

%% Constants (kept)
ME   = 9.10938356e-31;
MI   = 1.6737236e-27;
Q    = 1.60217662e-19;
EPS0 = 8.854187e-12;
amu  = 18;

%% Files
outnc  = 'profilesDIIID_196154_multi.nc';   % contains x,z, ni_all,u*_all, Z_all,q_all
ncfile = 'profilesDIIID_196154.nc';                % single-fluid (ne, te, ti, vr, vt, vz)

%% --- Read grids
R = ncread(ncfile,'x');         % [nR]
z = ncread(ncfile,'z');         % [nZ]
Rvec = R(:); zvec = z(:);
[Rvec, iR] = sort(Rvec,'ascend');
[zvec, iz] = sort(zvec,'ascend');

%% --- Read single-fluid fields
ne = ncread(ncfile,'ne');
te = ncread(ncfile,'te');      % electron temp (not used for carbon Ti)
vr = ncread(ncfile,'vr');
vt = ncread(ncfile,'vt');
vz = ncread(ncfile,'vz');

% deuterium / single-fluid Ti (this is what we will use for "Ti" everywhere)
ti = ncread(ncfile,'ti');

% sort onto (Rvec,zvec)
ne = ne(iR,iz); te = te(iR,iz);
vr = vr(iR,iz); vt = vt(iR,iz); vz = vz(iR,iz);
ti = ti(iR,iz);     % <- this is the Ti we'll use (from deuterium)
% if your single-fluid arrays were stored as (z,R), swap once:
if size(ne,1)==numel(zvec) && size(ne,2)==numel(Rvec)
    ne = ne.'; te = te.'; vr = vr.'; vt = vt.'; vz = vz.'; ti = ti.';
end

% %% --- Read multispecies metadata & fields
% Z_all = ncread(ncfile,'atomic_number');  % [ns]
% q_all = ncread(ncfile,'charge_number');  %#ok<NASGU>
% 
% ni_all_raw = ncread(ncfile,'ni_all');    % may not be [nR x nZ x ns] as-is
% uR_all_raw = ncread(ncfile,'uR_all');
% uZ_all_raw = ncread(ncfile,'uZ_all');
% uT_all_raw = ncread(ncfile,'uT_all');

%% ------------------- Read-back sanity plots -------------------
Rr = ncread(outnc,'x'); Zz = ncread(outnc,'z');
ne_m = ncread(outnc,'ne');
te_m = ncread(outnc,'te');
ni_all = ncread(outnc,'ni_all');

figure; imagesc(Zz, Rr, ne_m); set(gca,'YDir','normal','FontSize',12);
xlabel('Z [m]'); ylabel('R [m]'); colorbar; title('Masked n_e on (R,Z) from NetCDF');

figure; imagesc(Rr, Zz, ni_all(:,:,6)'); set(gca,'YDir','normal','FontSize',12);
xlabel('Z [m]'); ylabel('R [m]'); colorbar; title('Masked n_e on (R,Z) from NetCDF');
 hold on;
    plot(g.lim(1,:), g.lim(2,:), 'k', 'LineWidth', 1);
    contour(X, Y, psiN, [psiN_mask psiN_mask], 'k--', 'LineWidth', 1.2);

    % 2) Example species slice (e.g., k = 6) — no transpose needed
k = 6;
figure; imagesc(Zz, Rr, ni_all(:,:,k)); axis xy; set(gca,'FontSize',12);
xlabel('Z [m]'); ylabel('R [m]'); colorbar;
title(sprintf('Masked n_i (species %d) on (R,Z) from NetCDF', k));

hold on;
plot(g.lim(2,:), g.lim(1,:), 'k', 'LineWidth', 1.0);
contour(Zz, Rr, psi_m, [psiN_mask psiN_mask], 'k--', 'LineWidth', 1.2);
hold off;


%% --- Select ONLY Carbon (Z=6); sum all carbon charge states if multiple
idxC = find(Z_all==6);
if isempty(idxC)
    error('No carbon species (Z==6) found.');
end

% carbon density (sum all carbon charges)
niC = ni_all(:,:,9);

% carbon velocities (density-weighted over carbon charges)
tiny = 1e-60;
wC   = ni_all(:,:,idxC);                         % [nR x nZ x nC]
denC = sum(wC,3,'omitnan'); denC(denC<tiny)=tiny;
uR_C = sum(wC .* uR_all(:,:,idxC), 3, 'omitnan') ./ denC;
uZ_C = sum(wC .* uZ_all(:,:,idxC), 3, 'omitnan') ./ denC;
uT_C = sum(wC .* uT_all(:,:,idxC), 3, 'omitnan') ./ denC;

% IMPORTANT: "only Ti from deuterium"
% -> do NOT use any species ti; use single-fluid ti field already read/sorted above.

%% --- Read-back sanity style quick plots (match your axes/transposes)
% (1) electron density (single-fluid)
figure; imagesc(zvec, Rvec, ne); axis xy; set(gca,'FontSize',12);
xlabel('Z [m]'); ylabel('R [m]'); colorbar; title('Single-fluid n_e on (R,Z)');

% (2) carbon density niC: two equivalent ways (pick the one you prefer)
%     a) imagesc with x=R, y=Z by transposing
figure; imagesc(Rvec, zvec, niC.'); set(gca,'YDir','normal','FontSize',12);
xlabel('R [m]'); ylabel('Z [m]'); colorbar; title('Carbon (Z=6) n_i on (R,Z) [transpose method]');
hold on; plot(g.lim(1,:), g.lim(2,:), 'k', 'LineWidth', 1);

%     b) imagesc with axis xy, no transpose (x=z, y=R)
figure; imagesc(zvec, Rvec, niC); axis xy; set(gca,'FontSize',12);
xlabel('Z [m]'); ylabel('R [m]'); colorbar; title('Carbon (Z=6) n_i on (R,Z) [axis xy method]');
hold on; plot(g.lim(2,:), g.lim(1,:), 'k', 'LineWidth', 1);

%% --- Interpolants on Euclidean (R,Z)
interp_mode = 'linear';
extrap_mode = 'none';

F_niC = griddedInterpolant({Rvec,zvec}, niC, interp_mode, extrap_mode);
F_tiD = griddedInterpolant({Rvec,zvec}, ti,  interp_mode, extrap_mode);   % deuterium Ti
F_uR  = griddedInterpolant({Rvec,zvec}, uR_C, interp_mode, extrap_mode);
F_uZ  = griddedInterpolant({Rvec,zvec}, uZ_C, interp_mode, extrap_mode);
F_uT  = griddedInterpolant({Rvec,zvec}, uT_C, interp_mode, extrap_mode);

%% --- Map to your limiter geometry
r_centroid = sqrt(centroid(:,1).^2 + centroid(:,2).^2);
Rq = r_centroid;
Zq = centroid(:,3);

niC_surf = F_niC(Rq, Zq);           % carbon density on surface
tiD_surf = F_tiD(Rq, Zq);           % *deuterium* Ti on surface
uR_surf  = F_uR(Rq, Zq);
uZ_surf  = F_uZ(Rq, Zq);
uT_surf  = F_uT(Rq, Zq);
vC_surf  = sqrt(uR_surf.^2 + uZ_surf.^2 + uT_surf.^2);

%% --- 3D surface plots: Carbon n_i, Deuterium T_i, Carbon |u|
figure;
patch(transpose(X), transpose(Y), transpose(Z), niC_surf, ...
      'FaceAlpha', 1, 'EdgeAlpha', 0.25);
title('Carbon (Z=6) density on geometry');
xlabel('x [mm]'); ylabel('y [mm]'); zlabel('z [mm]');
axis equal; view(3); colorbar; colormap(parula);

figure;
patch(transpose(X), transpose(Y), transpose(Z), tiD_surf, ...
      'FaceAlpha', 1, 'EdgeAlpha', 0.25);
title('Deuterium T_i on geometry (used for carbon too)');
xlabel('x [mm]'); ylabel('y [mm]'); zlabel('z [mm]');
axis equal; view(3); colorbar; colormap(parula);

figure;
patch(transpose(X), transpose(Y), transpose(Z), vC_surf, ...
      'FaceAlpha', 1, 'EdgeAlpha', 0.25);
title('Carbon (Z=6) speed |u| on geometry');
xlabel('x [mm]'); ylabel('y [mm]'); zlabel('z [mm]');
axis equal; view(3); colorbar; colormap(parula);

%% Diagnostics
fprintf('Carbon(Z=6) |n_i| grid range: [%.3e, %.3e]\n', min(niC(:)), max(niC(:)));
fprintf('Carbon(Z=6) |n_i| surface  : [%.3e, %.3e]\n', min(niC_surf), max(niC_surf));

%% -------- helper to coerce species dims to [nR x nZ x nS]
function A = normalize_species_dims(Ain, nR, nZ, nS)
    sz = size(Ain);
    if numel(sz) ~= 3
        error('Species array must be 3D, got %s', mat2str(sz));
    end
    % If already correct
    if isequal(sz, [nR, nZ, nS]), A = Ain; return; end
    % Try candidate permutations
    cands = {
        [1 2 3], [nR nZ nS];
        [2 3 1], [nZ nS nR];
        [3 1 2], [nS nR nZ];
        [2 1 3], [nZ nR nS];
        [3 2 1], [nS nZ nR];
        [1 3 2], [nR nS nZ];
    };
    for i=1:size(cands,1)
        perm = cands{i,1}; tgt = cands{i,2};
        if all(sz(perm)==tgt)
            A = permute(Ain, perm);
            if isequal(size(A), [nR nZ nS]), return; end
        end
    end
    % Fallback: pick dims closest to nR/nZ, rest is nS
    [~,iR] = min(abs(sz - nR));
    [~,iZ] = min(abs(sz - nZ));
    iS = setdiff(1:3, [iR iZ], 'stable');
    A = permute(Ain, [iR iZ iS]);
    if ~isequal(size(A), [nR nZ nS])
        error('Could not coerce to [nR x nZ x nS], got %s', mat2str(size(A)));
    end
end