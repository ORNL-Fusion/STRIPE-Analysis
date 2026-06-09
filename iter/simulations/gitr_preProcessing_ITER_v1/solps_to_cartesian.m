%% SOLPS-ITER to Cartesian grid (D + Ne species only)
% Converts cell-centered SOLPS data (Geo.pr/pz, State.na, State.ne)
% to a uniform (R,Z) grid for visualization or coupling with GITR.

clear; clc; close all;

%% === INPUT ===
matFile = 'solps_iter.mat';           % Input file from SOLPS
outName = 'solps_iter_cartesian.mat'; % Output file

%% === LOAD DATA ===
S = load(matFile);
Geo   = S.Geo;
State = S.State;

% --- Geometry ---
r = double(Geo.pr);   % (4 x Nc)
z = double(Geo.pz);   % (4 x Nc)
Nc = size(r,2);

% --- Plasma state ---
ne_cv = double(State.ne(:));    % electron density (Nc x 1)
na_cv = double(State.na);       % all species densities (Nc x Ns)
ua_cv = [];
if isfield(State,'ua')
    ua_cv = double(State.ua);   % parallel flow (optional)
end

zn = double(State.zn(:));       % nuclear charge (Z)
am = double(State.am(:));       % atomic mass

fprintf('Loaded %d cells and %d species.\n', Nc, size(na_cv,2));

%% === IDENTIFY SPECIES (D + Ne) ===
idxD  = find(zn==1 & am==2);          % D0, D+
idxNe = find(zn==10 & am==20);        % Ne0..Ne10

idxD0    = idxD(1);
idxD1    = idxD(2);
idxNe0   = idxNe(1);
idxNeIon = idxNe(2:end);

fprintf('Detected %d neon charge states (Ne1+ to Ne%d+).\n', numel(idxNeIon), numel(idxNeIon));

%% === DEFINE CARTESIAN GRID ===
rmin = min(r(:)); rmax = max(r(:));
zmin = min(z(:)); zmax = max(z(:));
padR = 0.01*(rmax - rmin);
padZ = 0.01*(zmax - zmin);

rgrid = linspace(rmin+padR, rmax-padR, 400);
zgrid = linspace(zmin+padZ, zmax-padZ, 800);
[Rm, Zm] = meshgrid(rgrid, zgrid);

%% === HELPER FUNCTION (defined at bottom) ===
fill_by_cells = @(vals) paintCells(vals, r, z, Rm, Zm);

%% === ELECTRONS ===
disp('Mapping electron density...');
ne_q = fill_by_cells(ne_cv);

%% === DEUTERIUM ===
disp('Mapping D0 and D+ ...');
nD0_q = fill_by_cells(na_cv(:, idxD0));
nD1_q = fill_by_cells(na_cv(:, idxD1));

%% === NEON CHARGE STATES ===
disp('Mapping Ne species (Ne0..Ne10+) ...');
nNe0_q = fill_by_cells(na_cv(:, idxNe0));

nNe_q = cell(numel(idxNeIon),1);
for k = 1:numel(idxNeIon)
    nNe_q{k} = fill_by_cells(na_cv(:, idxNeIon(k)));
end

%% === TOTALS (robust) ===
disp('Computing totals ...');
% Force row vectors for safe concatenation
idxD1     = idxD1(:)'; 
idxNeIon  = idxNeIon(:)';

% Charged species = D+ + all Ne charge states
chargedIdx = [idxD1, idxNeIon];
chargedIdx = chargedIdx(chargedIdx <= size(na_cv,2));

% Remove invalid or NaN columns
validCols = all(isfinite(na_cv(:,chargedIdx)),1);
chargedIdx = chargedIdx(validCols);

% Sum charged ions
nIon_total_cv = sum(na_cv(:,chargedIdx), 2, 'omitnan');
nIon_total_q  = fill_by_cells(nIon_total_cv);

% Neon total (neutral + ions)
idxNe_all = [idxNe0, idxNeIon(:)'];
idxNe_all = idxNe_all(idxNe_all <= size(na_cv,2));
nNe_total_cv = sum(na_cv(:,idxNe_all), 2, 'omitnan');
nNe_total_q  = fill_by_cells(nNe_total_cv);

%% === PLOTS ===
figure; imagesc(rgrid,zgrid,ne_q);
set(gca,'YDir','normal','ColorScale','log');
colorbar; title('n_e [m^{-3}]'); xlabel('R [m]'); ylabel('Z [m]');
axis equal tight;

figure; imagesc(rgrid,zgrid,nD1_q);
set(gca,'YDir','normal','ColorScale','log');
colorbar; title('n_{D^+} [m^{-3}]'); xlabel('R [m]'); ylabel('Z [m]');
axis equal tight;

figure; imagesc(rgrid,zgrid,nNe0_q);
set(gca,'YDir','normal','ColorScale','log');
colorbar; title('n_{Ne^0} [m^{-3}]'); xlabel('R [m]'); ylabel('Z [m]');
axis equal tight;

if ~isempty(nNe_q)
    figure; imagesc(rgrid,zgrid,nNe_q{1});
    set(gca,'YDir','normal','ColorScale','log');
    colorbar; title('n_{Ne^{+}} [m^{-3}]'); xlabel('R [m]'); ylabel('Z [m]');
    axis equal tight;
end

%% === SAVE OUTPUT ===
save(outName, 'rgrid','zgrid','ne_q','nD0_q','nD1_q',...
    'nNe0_q','nNe_q','nNe_total_q','nIon_total_q','-v7.3');

fprintf('\n✅ Saved gridded output to: %s\n', outName);

%% === LOCAL FUNCTION ===
function grid_vals = paintCells(cv_vals, r4xN, z4xN, Rm, Zm)
    % Paints cell-centered values onto Cartesian mesh
    grid_vals = zeros(size(Rm));
    N = size(r4xN,2);
    for i = 1:N
        [in,on] = inpolygon(Rm, Zm, r4xN(:,i), z4xN(:,i));
        if any(in(:)) || any(on(:))
            grid_vals(in|on) = cv_vals(i);
        end
    end
end