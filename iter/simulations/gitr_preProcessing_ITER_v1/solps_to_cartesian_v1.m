%% SOLPS-ITER to Cartesian grid (D + Ne case, full profiles)
% Converts cell-centered SOLPS data (Geo.pr/pz, State.ne, State.na, State.te, State.ti, State.ua)
% into a uniform (R,Z) Cartesian grid with individual and summed quantities.

clear; clc; close all;

%% === INPUT ===
matFile = 'solps_iter.mat';
outName = 'solps_iter_cartesian_full.mat';

%% === LOAD DATA ===
S = load(matFile);
Geo   = S.Geo;
State = S.State;

% --- Geometry ---
r = double(Geo.pr);   % (4 x Nc)
z = double(Geo.pz);   % (4 x Nc)
Nc = size(r,2);

% --- Plasma state ---
ne_cv = double(State.ne(:));    % electron density
na_cv = double(State.na);       % species densities (Nc x Ns)
te_cv = double(State.te(:));    % electron temperature
ti_cv = double(State.ti(:));    % ion temperature
ua_cv = double(State.ua);       % parallel flows (Nc x Ns)

zn = double(State.zn(:));       % charge number
am = double(State.am(:));       % atomic mass
nSpecies = size(na_cv,2);

fprintf('Loaded %d cells, %d species, %d parallel flows.\n', Nc, nSpecies, size(ua_cv,2));

%% === IDENTIFY SPECIES ===
idxD  = find(zn==1 & am==2);     % D0, D+
idxNe = find(zn==10 & am==20);   % Ne0..Ne10

idxD0    = idxD(1);
idxD1    = idxD(2);
idxNe0   = idxNe(1);
idxNeIon = idxNe(2:end);

fprintf('Detected %d neon charge states (Ne1+ to Ne%d+)\n', numel(idxNeIon), numel(idxNeIon));

%% === DEFINE CARTESIAN GRID ===
rmin = min(r(:)); rmax = max(r(:));
zmin = min(z(:)); zmax = max(z(:));
padR = 0.01*(rmax - rmin);
padZ = 0.01*(zmax - zmin);

rgrid = linspace(rmin+padR, rmax-padR, 400);
zgrid = linspace(zmin+padZ, zmax-padZ, 800);
[Rm, Zm] = meshgrid(rgrid, zgrid);

%% === HELPER FUNCTION ===
fill_by_cells = @(vals) paintCells(vals, r, z, Rm, Zm);

%% === ELECTRONS ===
disp('Mapping n_e and T_e ...');
ne_q = fill_by_cells(ne_cv);
Te_q = fill_by_cells(te_cv);

%% === ION TEMPERATURE ===
disp('Mapping T_i ...');
Ti_q = fill_by_cells(ti_cv);

%% === DEUTERIUM DENSITIES ===
disp('Mapping D0 and D+ ...');
nD0_q = fill_by_cells(na_cv(:,idxD0));
nD1_q = fill_by_cells(na_cv(:,idxD1));

%% === NEON DENSITIES ===
disp('Mapping Ne0..Ne10+ ...');
nNe0_q = fill_by_cells(na_cv(:,idxNe0));
nNe_q = cell(numel(idxNeIon),1);
for k = 1:numel(idxNeIon)
    nNe_q{k} = fill_by_cells(na_cv(:,idxNeIon(k)));
end

%% === PARALLEL FLOWS (Ua) ===
disp('Mapping all 16 Ua parallel flows ...');
Ua_q = cell(1, nSpecies);
for s = 1:nSpecies
    Ua_q{s} = fill_by_cells(ua_cv(:,s));
end

%% === TOTALS ===
disp('Computing total charged and neon densities ...');

idxD1    = idxD1(:)'; 
idxNeIon = idxNeIon(:)';

chargedIdx = [idxD1, idxNeIon];
chargedIdx = chargedIdx(chargedIdx <= size(na_cv,2));
validCols = all(isfinite(na_cv(:,chargedIdx)),1);
chargedIdx = chargedIdx(validCols);

nIon_total_cv = sum(na_cv(:,chargedIdx), 2, 'omitnan');
nIon_total_q  = fill_by_cells(nIon_total_cv);

idxNe_all = [idxNe0, idxNeIon(:)'];
idxNe_all = idxNe_all(idxNe_all <= size(na_cv,2));
nNe_total_cv = sum(na_cv(:,idxNe_all), 2, 'omitnan');
nNe_total_q  = fill_by_cells(nNe_total_cv);

%% === PLOTS ===
figure('Name','Electron Density'); imagesc(rgrid,zgrid,ne_q);
set(gca,'YDir','normal','ColorScale','log'); colorbar;
title('n_e [m^{-3}]'); xlabel('R [m]'); ylabel('Z [m]'); axis equal tight;

figure('Name','Electron Temperature'); imagesc(rgrid,zgrid,Te_q);
set(gca,'YDir','normal','ColorScale','linear'); colorbar;
title('T_e [eV]'); xlabel('R [m]'); ylabel('Z [m]'); axis equal tight;

figure('Name','Ion Temperature'); imagesc(rgrid,zgrid,Ti_q);
set(gca,'YDir','normal','ColorScale','linear'); colorbar;
title('T_i [eV]'); xlabel('R [m]'); ylabel('Z [m]'); axis equal tight;

figure('Name','D+ Density'); imagesc(rgrid,zgrid,nD1_q);
set(gca,'YDir','normal','ColorScale','log'); colorbar;
title('n_{D^+} [m^{-3}]'); xlabel('R [m]'); ylabel('Z [m]'); axis equal tight;

figure('Name','Ne^0 Density'); imagesc(rgrid,zgrid,nNe0_q);
set(gca,'YDir','normal','ColorScale','log'); colorbar;
title('n_{Ne^0} [m^{-3}]'); xlabel('R [m]'); ylabel('Z [m]'); axis equal tight;

% Plot all Neon charge states
for k = 1:numel(nNe_q)
    figure('Name',sprintf('Ne^{%d+}',k));
    imagesc(rgrid,zgrid,nNe_q{k});
    set(gca,'YDir','normal','ColorScale','log');
    colorbar;
    title(sprintf('n_{Ne^{%d+}} [m^{-3}]',k));
    xlabel('R [m]'); ylabel('Z [m]');
    axis equal tight;
end

% Plot all 16 Ua flows
for s = 1:nSpecies
    figure('Name',sprintf('Ua Species %d',s));
    imagesc(rgrid,zgrid,Ua_q{s});
    set(gca,'YDir','normal','ColorScale','linear');
    colorbar;
    title(sprintf('U_{||, species %d} [m/s]',s));
    xlabel('R [m]'); ylabel('Z [m]');
    axis equal tight;
end

% Plot totals
figure('Name','Total Ion Density'); imagesc(rgrid,zgrid,nIon_total_q);
set(gca,'YDir','normal','ColorScale','log'); colorbar;
title('Total Charged Ion Density [m^{-3}]');
xlabel('R [m]'); ylabel('Z [m]'); axis equal tight;

figure('Name','Total Neon Density'); imagesc(rgrid,zgrid,nNe_total_q);
set(gca,'YDir','normal','ColorScale','log'); colorbar;
title('Total Neon Density [m^{-3}]');
xlabel('R [m]'); ylabel('Z [m]'); axis equal tight;

%% === SAVE ===
save(outName,'rgrid','zgrid','ne_q','Te_q','Ti_q',...
    'nD0_q','nD1_q','nNe0_q','nNe_q',...
    'Ua_q','nIon_total_q','nNe_total_q','-v7.3');

fprintf('\n✅ Saved full profiles and totals to: %s\n', outName);

%% === LOCAL FUNCTION ===
function grid_vals = paintCells(cv_vals, r4xN, z4xN, Rm, Zm)
    % Paints each cell's scalar value onto Cartesian mesh
    grid_vals = zeros(size(Rm));
    N = size(r4xN,2);
    for i = 1:N
        [in,on] = inpolygon(Rm, Zm, r4xN(:,i), z4xN(:,i));
        if any(in(:)) || any(on(:))
            grid_vals(in|on) = cv_vals(i);
        end
    end
end