%% Build extrapolated profiles with SOL-only continuation (200882)
close all; clear; clc;

%% --- Load SOLPS on native rS–zS grid
fileSOLPS = 'interpolated_values_196154.nc';
rS = ncread(fileSOLPS,'gridr');           % [Nr]
zS = ncread(fileSOLPS,'gridz');           % [Nz]
neS = ncread(fileSOLPS,'ne');             % [Nr x Nz]
TeS = ncread(fileSOLPS,'te');             % [Nr x Nz]
gradTiS = -ncread(fileSOLPS,'gradTi');
gradTirS =  ncread(fileSOLPS,'gradTir');  %#ok<NASGU> (kept for completeness)
vrS = ncread(fileSOLPS,'vr');
vtS = ncread(fileSOLPS,'vt');
vzS = ncread(fileSOLPS,'vz');

%% --- Floors / cleaning
n_min = 1e10;
T_min = 10;
neS(neS<=0 | ~isfinite(neS)) = NaN;
TeS(TeS<=0 | ~isfinite(TeS)) = NaN;
vrS(~isfinite(vrS)) = NaN; vzS(~isfinite(vzS)) = NaN; vtS(~isfinite(vtS)) = NaN;
gradTiS(~isfinite(gradTiS)) = NaN;
TeS(TeS < T_min) = T_min;

%% --- EFIT
read_efit_data;     % must populate struct g (limiter, axis, psirz info, etc.)

%% --- Quick looks on native grid
figure; imagesc(rS,zS,neS'); set(gca,'YDir','normal'); colorbar;
title('n_e (SOLPS native)'); hold on; plot(g.lim(1,:),g.lim(2,:),'r');

figure; imagesc(rS,zS,TeS'); set(gca,'YDir','normal'); colorbar;
title('T_e (SOLPS native)'); hold on; plot(g.lim(1,:),g.lim(2,:),'r');

%% --- Build X–Y grid bounded by limiter
num_points = 1000;
[X,Y] = meshgrid( linspace(min(g.lim(1,:)),max(g.lim(1,:)),num_points), ...
                  linspace(min(g.lim(2,:)),max(g.lim(2,:)),num_points) );

%% --- \psi_N on X–Y (mask will be based on this)
[psiN_vec,~] = calc_psiN(g, X(:).', Y(:).', 0);
psiN_XY = reshape(psiN_vec, size(X));
mask_SOL = psiN_XY > 1.0;                         % only extrapolate here
inside_limiter = inpolygon(X, Y, g.lim(1,:), g.lim(2,:));

%% --- Interpolate original SOLPS to X–Y (kept inside \psi_N<=1)
base_ne = interp2(rS, zS, neS', X, Y, 'linear', NaN);
base_Te = interp2(rS, zS, TeS', X, Y, 'linear', NaN);
base_ne(~isfinite(base_ne)) = 0;
base_Te(~isfinite(base_Te)) = T_min;

%% --- Prepare native sample cloud (for occasional use)
[ZZ,RR] = meshgrid(zS,rS);         % [Nz x Nr]
rkron = kron(rS, ones(size(zS))')';
zkron = kron(zS, ones(size(rS)));
idx_valid = find(~isnan(neS'));
rcoords = rkron(idx_valid); zcoords = zkron(idx_valid);
coords  = [rcoords, zcoords];             % M×2 (r,z)
ok_ne   = neS'; ok_ne = ok_ne(idx_valid);
ok_Te   = TeS'; ok_Te = ok_Te(idx_valid);

%% --- Estimate decay slope near LCFS from midplane data (automatic)
R_line = rS(:).';  Z_line = zeros(size(R_line));
[psi_mid,~] = calc_psiN(g, R_line, Z_line, 0);
ne_mid = interp2(rS, zS, neS', R_line, Z_line, 'linear', NaN);
Te_mid = interp2(rS, zS, TeS', R_line, Z_line, 'linear', NaN);

band = psi_mid>0.95 & psi_mid<1.05 & ne_mid>0 & isfinite(ne_mid);
if nnz(band)>=5, p_ne = polyfit(R_line(band), log(ne_mid(band)), 1);
else,            p_ne = polyfit(R_line(isfinite(ne_mid)), log(ne_mid(isfinite(ne_mid))), 1);
end
bandT = psi_mid>0.95 & psi_mid<1.05 & Te_mid>0 & isfinite(Te_mid);
if nnz(bandT)>=5, p_Te = polyfit(R_line(bandT), log(Te_mid(bandT)), 1);
else,             p_Te = polyfit(R_line(isfinite(Te_mid)), log(Te_mid(isfinite(Te_mid))), 1);
end
Ln  = max(0.5*abs(1/p_ne(1)), 0.005);   % clamp to sane mins
LTe = max(0.5*abs(1/p_Te(1)), 0.005);

%% --- Extract LCFS contour (psiN = 1) on X–Y
C = contourc(X(1,:), Y(:,1), psiN_XY, [1 1]);
k = 1; nseg = C(2,1);
lcfs_xy = C(:, 2:1+nseg);
R_lcfs = lcfs_xy(1,:).'; Z_lcfs = lcfs_xy(2,:).';

%% --- Edge values just inside LCFS (to anchor extrapolation)
edge_band = (psiN_XY >= 0.98) & (psiN_XY <= 1.02);
edge_band(psiN_XY > 1) = false;      % only inside
if nnz(edge_band) < 100
    edge_band = psiN_XY <= 1.0;      % fallback
end
edge_R = X(edge_band); edge_Z = Y(edge_band);
edge_ne_vals = base_ne(edge_band);
edge_Te_vals = base_Te(edge_band);

F_ne_edge = scatteredInterpolant(edge_R, edge_Z, edge_ne_vals, 'natural', 'nearest');
F_Te_edge = scatteredInterpolant(edge_R, edge_Z, edge_Te_vals, 'natural', 'nearest');

%% --- Distance to LCFS (approx normal distance)
P = [X(:) Y(:)];
Q = [R_lcfs Z_lcfs];
idxN = knnsearch(Q, P);                    % nearest LCFS sample
d_RZ = vecnorm(P - Q(idxN,:), 2, 2);       % unsigned distance
d_RZ = reshape(d_RZ, size(X));
d_signed = d_RZ;                           % sign by psiN
d_signed(psiN_XY < 1) = -d_signed(psiN_XY < 1);

%% --- Extrapolated fields (outside only), edge-anchored exponential
ne_edge = F_ne_edge(X, Y);
Te_edge = F_Te_edge(X, Y);
ne_extrap = ne_edge .* exp( -max(d_signed,0) / Ln );
Te_extrap = max( T_min, Te_edge .* exp( -max(d_signed,0) / LTe ) );

%% --- Smooth blend across LCFS to avoid a hard wall
Delta = 0.02;                                 % width in psiN
alpha = min(max((psiN_XY - 1)/Delta, 0), 1);  % 0 inside, 1 outside

val_ne = (1-alpha).*base_ne + alpha.*ne_extrap;
val_Te = (1-alpha).*base_Te + alpha.*Te_extrap;

% enforce limiter
val_ne(~inside_limiter) = 0;
val_Te(~inside_limiter) = T_min;

%% --- Visual checks
figure; imagesc([min(X(:)) max(X(:))],[min(Y(:)) max(Y(:))], psiN_XY);
set(gca,'YDir','normal'); colorbar; title('\psi_N on X–Y'); hold on; plot(g.lim(1,:),g.lim(2,:),'k');

figure; imagesc([min(X(:)) max(X(:))],[min(Y(:)) max(Y(:))], val_ne);
set(gca,'YDir','normal'); colorbar; set(gca,'ColorScale','log');
title('n_e with SOL-only extrapolation'); hold on; plot(g.lim(1,:),g.lim(2,:),'r');

figure; imagesc([min(X(:)) max(X(:))],[min(Y(:)) max(Y(:))], val_Te);
set(gca,'YDir','normal'); colorbar; set(gca,'ColorScale','log');
title('T_e with SOL-only extrapolation'); hold on; plot(g.lim(1,:),g.lim(2,:),'r');

%% --- Your original velocity/gradTi extrapolation on native grid
inside_lim_native = inpolygon(RR, ZZ, g.lim(1,:), g.lim(2,:));
plasma_mask = neS > n_min;
decay_length = 0.2;

for iz = 1:length(zS)
    edge_idx = find(plasma_mask(:,iz), 1, 'last');
    if ~isempty(edge_idx) && edge_idx < length(rS)
        R_edge = rS(edge_idx);
        evr = vrS(edge_idx,iz); evz = vzS(edge_idx,iz);
        evt = vtS(edge_idx,iz); egT = gradTiS(edge_idx,iz);
        for ir = edge_idx+1:length(rS)
            if inside_lim_native(iz,ir)
                fac = exp(-(rS(ir)-R_edge)/decay_length);
                vrS(ir,iz)     = evr * fac;
                vzS(ir,iz)     = evz * fac;
                vtS(ir,iz)     = evt * fac;
                gradTiS(ir,iz) = egT * fac;
            end
        end
    end
end

%% --- Interpolate those to X–Y
val_vr     = interp2(rS, zS, vrS',     X, Y, 'linear');
val_vz     = interp2(rS, zS, vzS',     X, Y, 'linear');
val_vt     = interp2(rS, zS, vtS',     X, Y, 'linear');
val_gradTi = interp2(rS, zS, gradTiS', X, Y, 'linear');

%% --- Multi-field quick views
vars2D    = {val_ne, val_Te, val_vr, val_vz, val_vt, val_gradTi};
varNames2D= {'Electron Density (n_e)','Electron Temperature (T_e)', ...
             'Radial Velocity (v_r)','Poloidal Velocity (v_z)', ...
             'Toroidal Velocity (v_t)','Ion Temperature Gradient (\nablaT_i)'};

for i = 1:numel(vars2D)
    figure;
    imagesc([min(X(:)) max(X(:))],[min(Y(:)) max(Y(:))],vars2D{i});
    set(gca,'YDir','normal'); colorbar; set(gca,'ColorScale','linear');
    title(['Extrapolated ',varNames2D{i}]); hold on; plot(g.lim(1,:),g.lim(2,:),'r');
end

%% --- Save
save('extrapolated_data_200882_maskedSOL.mat', ...
     'X','Y','val_ne','val_Te','val_vr','val_vz','val_vt','val_gradTi', ...
     'psiN_XY','g','rS','zS');

disp('Done: kept SOLPS inside (psiN<=1), extrapolated only for psiN>1 with LCFS distance + smooth blend.');