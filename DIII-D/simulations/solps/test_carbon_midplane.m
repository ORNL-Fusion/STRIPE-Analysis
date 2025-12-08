%% Carbon (or any species) midplane plot — same method as deuterium

% --- Inputs needed from your script context ---
% fileSOLPS, rS [nr], zS [nz], T_min, p_ne, p_Te (from electron midplane fits)
% multispecies fields read as in your run:
%   ni_all: [ns x nr x nz] or [nz x nr x ns] (we'll permute)
%   ti_all: [ns x nr x nz] or [nz x nr x ns]
%   atomic_number: Z_all [ns], charge_number: q_all [ns]
close all; clear all;
fileSOLPS = 'interpolated_values.nc';
rS = ncread(fileSOLPS, 'gridr');    % Radial grid (Nr)
zS = ncread(fileSOLPS, 'gridz');    % Vertical grid (Nz)
neS = ncread(fileSOLPS, 'ne');      % Electron density (Nr x Nz)
TeS = ncread(fileSOLPS, 'te');      % Electron temperature (Nr x Nz)
gradTiS=-ncread(fileSOLPS,'gradTi');
gradTirS=ncread(fileSOLPS,'gradTir');
vrS = ncread(fileSOLPS,'vr');
vtS = ncread(fileSOLPS,'vt');
vzS = ncread(fileSOLPS,'vz');


Z_all = ncread(fileSOLPS,'atomic_number');
q_all = ncread(fileSOLPS,'charge_number');
ni_raw = ncread(fileSOLPS,'ni_all');
ti_raw = ncread(fileSOLPS,'ti_all');

nr = numel(rS); nz = numel(zS); ns = numel(Z_all);

% ---- read raw arrays from the netCDF ----
Z_all = ncread(fileSOLPS,'atomic_number');
q_all = ncread(fileSOLPS,'charge_number');
ni_raw = ncread(fileSOLPS,'ni_all');   % some ordering of ns,nr,nz
ti_raw = ncread(fileSOLPS,'ti_all');   % same ordering as ni_all

nr = numel(rS);
nz = numel(zS);
ns = numel(Z_all);

% ---- normalize to [ns x nr x nz] with a safe helper ----
ni_all = to_ns_nr_nz(ni_raw, ns, nr, nz);
ti_all = to_ns_nr_nz(ti_raw, ns, nr, nz);

% (then continue with your k_here selection and plotting...)

% ================= helper =================
function Aout = to_ns_nr_nz(Ain, ns, nr, nz)
    sz = size(Ain);
    if numel(sz) ~= 3
        error('Expected 3D array, got size %s', mat2str(sz));
    end

    % Exact matches of the 3 likely layouts we've seen:
    if isequal(sz, [ns, nr, nz])
        Aout = Ain;                       % already [ns x nr x nz]
    elseif isequal(sz, [nr, nz, ns])
        Aout = permute(Ain, [3, 1, 2]);   % [nr x nz x ns] -> [ns x nr x nz]
    elseif isequal(sz, [nz, nr, ns])
        Aout = permute(Ain, [3, 2, 1]);   % [nz x nr x ns] -> [ns x nr x nz]
    elseif isequal(sz, [ns, nz, nr])
        Aout = permute(Ain, [1, 3, 2]);   % [ns x nz x nr] -> [ns x nr x nz]
    elseif isequal(sz, [nr, ns, nz])
        Aout = permute(Ain, [2, 1, 3]);   % [nr x ns x nz] -> [ns x nr x nz]
    elseif isequal(sz, [nz, ns, nr])
        Aout = permute(Ain, [2, 3, 1]);   % [nz x ns x nr] -> [ns x nr x nz]
    else
        error('Dont know how to permute from size %s to [ns x nr x nz] with ns=%d, nr=%d, nz=%d.', ...
              mat2str(sz), ns, nr, nz);
    end
end

% ---- choose species: example carbon with any charge (pick one by index) ---
% If you know the species index you want, set k_here = <index>;
% otherwise, pick the first Z=6:
cand = find(Z_all==6);
if isempty(cand)
    error('No species with Z=6 found in atomic_number.');
end
k_here = cand(1);   % change if you want a different C charge state

% Extract that species as [nr x nz]
ni_k = squeeze(ni_all(k_here,:,:));   % [nr x nz]
Ti_k = squeeze(ti_all(k_here,:,:));   % [nr x nz]

% ---------- midplane extraction (robust) ----------
% Use the z grid location closest to 0; or optionally average a small band.
[~, iz0] = min(abs(zS - 0));
z_mid = zS(iz0);

% Midplane sampling radii (like your D block):
mpfx = linspace(2.14333, 2.20675, 1000); % fit window near separatrix
mpfy = mpfx*0 + z_mid;                   % use z closest to 0, not exactly 0

mpx  = linspace(2.14333, 2.5, 1000);     % full radial range for blended curve

% Interpolate along the midplane line:
ni_mid  = interp2(rS, zS, ni_k', mpfx, mpfy, 'linear', NaN);
Ti_mid  = interp2(rS, zS, Ti_k', mpfx, mpfy, 'linear', NaN);

% Keep only positive finite samples for fitting
good_n = isfinite(ni_mid) & (ni_mid > 0);
good_T = isfinite(Ti_mid) & (Ti_mid > 0);

% Fit slopes EXACTLY like D; fallback to electron slopes if too few points
if nnz(good_n) >= 3
    p_ni = polyfit(mpfx(good_n), log(ni_mid(good_n)), 1);
else
    warning('Too few midplane n_i points for species k=%d; using p_ne slope.', k_here);
    p_ni = p_ne;   % from your electron fit
end
if nnz(good_T) >= 3
    p_Ti = polyfit(mpfx(good_T), log(Ti_mid(good_T)), 1);
else
    warning('Too few midplane T_i points for species k=%d; using p_Te slope.', k_here);
    p_Ti = p_Te;   % from your electron fit
end

% Baseline profiles at midplane z=z_mid across full mpx (like D):
ni_base = interp2(rS, zS, ni_k', mpx, 0*mpx+z_mid, 'linear', NaN);
Ti_base = interp2(rS, zS, Ti_k', mpx, 0*mpx+z_mid, 'linear', NaN);
ni_base(~isfinite(ni_base)) = 0;
Ti_base(~isfinite(Ti_base)) = T_min;

% Same blending function as your D script:
interpfn = (mpx - 2.14333) / (2.20675 - 2.14333);
interpfn = min(max(interpfn, 0), 1);

% Extrapolated 1D curves (SAME formulae as D)
extrapolatedni1dC = interpfn .* exp(p_ni(2) + mpx * p_ni(1)) + (1 - interpfn) .* ni_base;
extrapolatedTi1dC = max(interpfn .* exp(p_Ti(2) + mpx * p_Ti(1)) + (1 - interpfn) .* Ti_base, T_min);

% For plotting the "midplane" raw curves over the same full mpx domain:
% (Use nearest neighbor to avoid NaNs if the midplane line leaves the mesh)
ni_mid_full = interp2(rS, zS, ni_k', mpx, 0*mpx+z_mid, 'nearest', NaN);
Ti_mid_full = interp2(rS, zS, Ti_k', mpx, 0*mpx+z_mid, 'nearest', NaN);
ni_mid_full(~isfinite(ni_mid_full)) = 0;
Ti_mid_full(~isfinite(Ti_mid_full)) = T_min;

% ------------- Plot (exact style you asked for) -------------
figure;
semilogy(mpx, ni_mid_full, 'b');  hold on;                        % "Density Midplane"
semilogy(mpx, extrapolatedni1dC, 'b.');                           % "Extrapolated Density"
semilogy(mpx, Ti_mid_full, 'r');                                  % "Temperature Midplane"
semilogy(mpx, extrapolatedTi1dC, 'r.');                           % "Extrapolated Temperature"
legend('Carbon Density Midplane', 'Extrapolated Carbon Density', ...
       'Carbon Temperature Midplane', 'Extrapolated Carbon Temperature', ...
       'Location','best');
xlabel('R [m]'); ylabel('Value'); set(gca,'YScale','log');
title(sprintf('Carbon Midplane Density and Temperature Extrapolation (species k=%d, Z=%d, q=%d)', ...
      k_here, Z_all(k_here), q_all(k_here)));
grid on;

