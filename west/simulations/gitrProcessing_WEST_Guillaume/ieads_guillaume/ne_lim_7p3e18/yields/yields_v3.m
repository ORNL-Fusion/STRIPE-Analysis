%% WEST surface erosion post-processing from ieads_comsol_* folders
% Includes:
% 1) pads missing yield / target entries with zeros
% 2) reads yields as yields(:,c), where c is charge state
% 3) multiplies ion density by species concentration
% 4) computes erosion flux, surface-integrated erosion rate, mass erosion rate
% 5) computes cumulative erosion CDF over surface faces

clear; clc; close all;

%% ======================= USER SETTINGS =======================
baseDir = '..';

speciesFolders = { ...
    'ieads_comsol_b1+', ...
    'ieads_comsol_b2+', ...
    'ieads_comsol_b3+', ...
    'ieads_comsol_c2+', ...
    'ieads_comsol_c3+', ...
    'ieads_comsol_c4+', ...
    'ieads_comsol_W2+'};

speciesLabels = { ...
    'B1+', 'B2+', 'B3+', ...
    'C2+', 'C3+', 'C4+', ...
    'W2+'};

speciesConc = [ ...
    0.12, ...
    0.12, ...
    0.16, ...
    0.06, ...
    0.06, ...
    0.08, ...
    0.001 ...
];

% Species/projectile atomic masses in amu
amu_species = [ ...
    10.81, ...    % B1+
    10.81, ...    % B2+
    10.81, ...    % B3+
    12.011, ...   % C2+
    12.011, ...   % C3+
    12.011, ...   % C4+
    183.84 ...    % W2+
];

amu_to_kg = 1.66053906660e-27;

density_col  = 11;
velocity_col = 5;

hasHeader = true;

if numel(speciesConc) ~= numel(speciesLabels)
    error('speciesConc must match speciesLabels.');
end

if numel(amu_species) ~= numel(speciesLabels)
    error('amu_species must match speciesLabels.');
end

%% ======================= FIND GEOMETRY =======================
geomFile = '';
for k = 1:numel(speciesFolders)
    cand = fullfile(baseDir, speciesFolders{k}, 'gitrGeometryPointPlane3d.cfg');
    if isfile(cand)
        geomFile = cand;
        break;
    end
end

if isempty(geomFile)
    error('Could not find gitrGeometryPointPlane3d.cfg.');
end

fprintf('Using geometry file:\n  %s\n\n', geomFile);

%% ======================= LOAD GEOMETRY =======================
fid = fopen(geomFile, 'r');
if fid == -1
    error('Failed to open geometry file: %s', geomFile);
end

lines = {};
while ~feof(fid)
    lines{end+1} = fgetl(fid);
end
fclose(fid);

for i = 3:numel(lines)
    tline = strtrim(lines{i});
    if isempty(tline)
        continue;
    end
    try
        evalc(tline);
    catch
    end
end

if ~exist('x1','var') || ~exist('x2','var') || ~exist('x3','var') || ...
   ~exist('y1','var') || ~exist('y2','var') || ~exist('y3','var') || ...
   ~exist('z1','var') || ~exist('z2','var') || ~exist('z3','var')
    error('Geometry variables x1..z3 were not created.');
end

subset = 1:length(x1);

X = [x1(subset).', x2(subset).', x3(subset).'];
Y = [y1(subset).', y2(subset).', y3(subset).'];
Z = [z1(subset).', z2(subset).', z3(subset).'];

Nfaces = size(X,1);
surf_inds = 1:Nfaces;

fprintf('Detected %d surface faces from geometry.\n\n', Nfaces);

%% ======================= FACE AREAS ==========================
p1 = [X(:,1), Y(:,1), Z(:,1)];
p2 = [X(:,2), Y(:,2), Z(:,2)];
p3 = [X(:,3), Y(:,3), Z(:,3)];

area = 0.5 .* vecnorm(cross(p2 - p1, p3 - p1, 2), 2, 2);
area(~isfinite(area)) = 0;

fprintf('Total surface area = %.6e m^2\n\n', sum(area));

%% ======================= INITIALIZE TOTALS ===================
Ns = numel(speciesFolders);

total_ero_flux  = zeros(Nfaces,1);
total_ion_flux  = zeros(Nfaces,1);
density_total   = zeros(Nfaces,1);
yield_total_sum = zeros(Nfaces,1);

total_erosion = zeros(Nfaces,1);

erosion_rate_species      = zeros(Ns,1);   % particles/s or atoms/s proxy
mass_erosion_rate_species = zeros(Ns,1);   % kg/s

all_yield       = zeros(Nfaces,Ns);
all_flux        = zeros(Nfaces,Ns);
all_nvf         = zeros(Nfaces,Ns);
all_density     = zeros(Nfaces,Ns);
all_density_raw = zeros(Nfaces,Ns);
all_erosion     = zeros(Nfaces,Ns);

%% ======================= SPECIES LOOP ========================
for i = 1:Ns
    folder = fullfile(baseDir, speciesFolders{i});
    label  = speciesLabels{i};
    conc_i = speciesConc(i);

    targetPattern = fullfile(folder, 'Targets_*.txt');
    yieldPattern  = fullfile(folder, 'yields_*.csv');

    tFiles = dir(targetPattern);
    yFiles = dir(yieldPattern);

    if isempty(tFiles)
        warning('Skipping %s: no Targets_*.txt found.', folder);
        continue;
    end

    if isempty(yFiles)
        warning('Skipping %s: no yields_*.csv found.', folder);
        continue;
    end

    targetFile = fullfile(folder, tFiles(1).name);
    yieldFile  = fullfile(folder, yFiles(1).name);

    fprintf('Reading %s\n', label);
    fprintf('  Targets      : %s\n', targetFile);
    fprintf('  Yields       : %s\n', yieldFile);
    fprintf('  Concentration: %.6g\n', conc_i);

    data = load_numeric_targets(targetFile, hasHeader);

    if size(data,2) < max(density_col, velocity_col)
        error('File %s has only %d columns.', targetFile, size(data,2));
    end

    charge_state = parse_charge_state(label);
    Y_eff = load_yield_vector(yieldFile, Nfaces, charge_state);

    dens_raw = pad_or_trim_vector(data(:, density_col), Nfaces, ...
        sprintf('Density column in %s', targetFile));

    dens_data = conc_i .* dens_raw;

    v_data = pad_or_trim_vector(data(:, velocity_col), Nfaces, ...
        sprintf('Velocity column in %s', targetFile));

    dens_data(~isfinite(dens_data)) = 0;
    v_data(~isfinite(v_data))       = 0;
    Y_eff(~isfinite(Y_eff))         = 0;

    ion_flux_species = abs(dens_data .* v_data);
    eroded_flux = abs(Y_eff .* ion_flux_species);

    writematrix(eroded_flux, sprintf('eroded_flux_%s.txt', label));

    erosion = eroded_flux .* area;
    erosion_surface = erosion;

    erosion_rate_species(i) = sum(erosion_surface, 'omitnan');
    mass_erosion_rate_species(i) = erosion_rate_species(i) .* ...
        amu_species(i) .* amu_to_kg;

    total_erosion = total_erosion + erosion(surf_inds);

    total_ion_flux  = total_ion_flux  + ion_flux_species;
    total_ero_flux  = total_ero_flux  + eroded_flux;
    density_total   = density_total   + dens_data;
    yield_total_sum = yield_total_sum + Y_eff;

    all_yield(:,i)       = Y_eff;
    all_flux(:,i)        = eroded_flux;
    all_nvf(:,i)         = ion_flux_species;
    all_density(:,i)     = dens_data;
    all_density_raw(:,i) = dens_raw;
    all_erosion(:,i)     = erosion;

    fprintf('  Erosion rate      = %.6e atoms/s\n', erosion_rate_species(i));
    fprintf('  Mass erosion rate = %.6e kg/s\n\n', mass_erosion_rate_species(i));

    figure('Name',['Species: ', label], 'Color','w');

    subplot(4,1,1);
    plot_linear_patch(X,Y,Z,Y_eff, ['Yield: ', label]);

    subplot(4,1,2);
    plot_log_patch(X,Y,Z,ion_flux_species, ['Ion flux: ', label]);

    subplot(4,1,3);
    plot_log_patch(X,Y,Z,eroded_flux, ['Eroded flux: ', label]);

    subplot(4,1,4);
    plot_log_patch(X,Y,Z,erosion, ['Surface erosion rate: ', label]);
end

%% ======================= TOTAL EROSION CDF ===================
erosion_inds = find(total_erosion);
erosion_sub = total_erosion(erosion_inds);

if isempty(erosion_sub)
    warning('total_erosion is zero everywhere. CDF cannot be computed.');
    erosion_sub_cdf = [];
    erosion_rate = 0;
else
    erosion_sub_cdf = cumsum(erosion_sub);
    erosion_rate = erosion_sub_cdf(end);
    erosion_sub_cdf = erosion_sub_cdf ./ erosion_sub_cdf(end);
end

mass_erosion_rate_total = sum(mass_erosion_rate_species, 'omitnan');

fprintf('\n===== TOTAL EROSION SUMMARY =====\n');
fprintf('Total erosion rate      = %.6e atoms/s\n', erosion_rate);
fprintf('Total mass erosion rate = %.6e kg/s\n', mass_erosion_rate_total);
fprintf('Total mass erosion rate = %.6e g/s\n', mass_erosion_rate_total * 1e3);

%% ======================= EFFECTIVE YIELDS ====================
yield_eff_total = zeros(Nfaces,1);
mask = total_ion_flux ~= 0;
yield_eff_total(mask) = total_ero_flux(mask) ./ total_ion_flux(mask);

yield_eff_B = zeros(Nfaces,1);
yield_eff_C = zeros(Nfaces,1);
yield_eff_W = zeros(Nfaces,1);

b_idx = find(startsWith(speciesLabels, 'B'));
c_idx = find(startsWith(speciesLabels, 'C'));
w_idx = find(startsWith(speciesLabels, 'W'));

if ~isempty(b_idx)
    num_B = sum(all_yield(:,b_idx) .* all_nvf(:,b_idx), 2, 'omitnan');
    den_B = sum(all_nvf(:,b_idx), 2, 'omitnan');
    mB = den_B ~= 0;
    yield_eff_B(mB) = num_B(mB) ./ den_B(mB);
end

if ~isempty(c_idx)
    num_C = sum(all_yield(:,c_idx) .* all_nvf(:,c_idx), 2, 'omitnan');
    den_C = sum(all_nvf(:,c_idx), 2, 'omitnan');
    mC = den_C ~= 0;
    yield_eff_C(mC) = num_C(mC) ./ den_C(mC);
end

if ~isempty(w_idx)
    num_W = sum(all_yield(:,w_idx) .* all_nvf(:,w_idx), 2, 'omitnan');
    den_W = sum(all_nvf(:,w_idx), 2, 'omitnan');
    mW = den_W ~= 0;
    yield_eff_W(mW) = num_W(mW) ./ den_W(mW);
end

%% ======================= TOTAL PLOTS =========================
figure('Name','WEST Total Erosion Flux','Color','w');

subplot(2,1,1);
plot_linear_patch(X,Y,Z,total_ero_flux, 'Total eroded flux');

subplot(2,1,2);
plot_log_patch(X,Y,Z,total_ero_flux, 'Total eroded flux log10');

figure('Name','WEST Total Surface Erosion Rate','Color','w');

subplot(2,1,1);
plot_linear_patch(X,Y,Z,total_erosion, 'Total surface erosion rate');

subplot(2,1,2);
plot_log_patch(X,Y,Z,total_erosion, 'Total surface erosion rate log10');

figure('Name','WEST Total Ion Flux and Density','Color','w');

subplot(2,1,1);
plot_log_patch(X,Y,Z,total_ion_flux, 'Total ion flux');

subplot(2,1,2);
plot_log_patch(X,Y,Z,density_total, 'Total density');

figure('Name','WEST Effective Yields','Color','w');

subplot(1,4,1);
plot_linear_patch(X,Y,Z,yield_eff_B, 'Effective yield - B');

subplot(1,4,2);
plot_linear_patch(X,Y,Z,yield_eff_C, 'Effective yield - C');

subplot(1,4,3);
plot_linear_patch(X,Y,Z,yield_eff_W, 'Effective yield - W');

subplot(1,4,4);
plot_linear_patch(X,Y,Z,yield_eff_total, 'Effective yield - Total');

if ~isempty(erosion_sub_cdf)
    figure('Name','Surface Erosion CDF','Color','w');
    plot(erosion_sub_cdf, 'LineWidth', 2);
    grid on;
    xlabel('Surface element index');
    ylabel('Cumulative erosion fraction');
    title('CDF of total surface erosion');
end

%% ======================= SAVE OUTPUTS ========================
writematrix(total_ero_flux,  'west_total_eroded_flux.txt');
writematrix(total_ion_flux,  'west_total_ion_flux.txt');
writematrix(total_erosion,   'west_total_surface_erosion_rate.txt');
writematrix(yield_eff_total, 'west_total_effective_yield.txt');
writematrix(yield_total_sum, 'west_total_yield_sum.txt');
writematrix(density_total,   'west_density_total_all_species.txt');

writematrix(all_yield,       'west_species_yields.txt');
writematrix(all_flux,        'west_species_eroded_flux.txt');
writematrix(all_nvf,         'west_species_ion_flux.txt');
writematrix(all_density,     'west_species_density_conc_weighted.txt');
writematrix(all_density_raw, 'west_species_density_raw.txt');
writematrix(all_erosion,     'west_species_surface_erosion_rate.txt');

writematrix(speciesConc(:),  'west_species_concentrations.txt');
writematrix(amu_species(:),  'west_species_amu.txt');

writematrix(erosion_rate_species,      'west_erosion_rate_species_atoms_per_s.txt');
writematrix(mass_erosion_rate_species, 'west_mass_erosion_rate_species_kg_per_s.txt');

writematrix(erosion_rate,              'west_total_erosion_rate_atoms_per_s.txt');
writematrix(mass_erosion_rate_total,   'west_total_mass_erosion_rate_kg_per_s.txt');

writematrix(erosion_inds,     'west_erosion_cdf_indices.txt');
writematrix(erosion_sub,      'west_erosion_cdf_values.txt');
writematrix(erosion_sub_cdf,  'west_erosion_cdf_normalized.txt');

fprintf('\nDone.\n');

%% ======================= LOCAL FUNCTIONS =====================
function vals = load_numeric_targets(fname, hasHeader)
    if hasHeader
        T = readtable(fname, 'Delimiter', ',', 'ReadVariableNames', true);
        vals = table2array(T);
    else
        vals = readmatrix(fname);
    end

    vals = double(vals);
    vals(~isfinite(vals)) = 0;
end

function c = parse_charge_state(label)
    tok = regexp(label, '([A-Za-z]+)(\d+)\+', 'tokens', 'once');

    if isempty(tok)
        error('Could not parse charge state from label: %s', label);
    end

    c = str2double(tok{2});

    if ~isfinite(c) || c < 1
        error('Invalid charge state parsed from label: %s', label);
    end
end

function y = load_yield_vector(fname, Nfaces, charge_state)
    yraw = readmatrix(fname);
    yraw = double(yraw);

    if isempty(yraw)
        warning('Yield file %s is empty. Filling with zero.', fname);
        y = zeros(Nfaces,1);
        return;
    end

    if isvector(yraw)
        yraw = yraw(:);
    end

    yraw(~isfinite(yraw)) = 0;

    if charge_state > size(yraw,2)
        warning('Yield file %s has only %d columns. Requested charge state %d. Using zeros.', ...
            fname, size(yraw,2), charge_state);
        ycol = zeros(0,1);
    else
        ycol = yraw(:, charge_state);
    end

    ycol = double(ycol(:));
    ycol(~isfinite(ycol)) = 0;

    if numel(ycol) == Nfaces + 1 && abs(ycol(1)) < 1e-14
        ycol = ycol(2:end);
    end

    y = zeros(Nfaces,1);
    ncopy = min(Nfaces, numel(ycol));
    y(1:ncopy) = ycol(1:ncopy);

    if numel(ycol) < Nfaces
        warning('Yield column c=%d in %s has %d values, expected %d. Padding zeros.', ...
            charge_state, fname, numel(ycol), Nfaces);
    elseif numel(ycol) > Nfaces
        warning('Yield column c=%d in %s has %d values, expected %d. Truncating.', ...
            charge_state, fname, numel(ycol), Nfaces);
    end

    y(~isfinite(y)) = 0;
end

function v = pad_or_trim_vector(v, Nfaces, nameStr)
    v = double(v(:));
    v(~isfinite(v)) = 0;

    if numel(v) > Nfaces
        warning('%s has %d entries, expected %d. Truncating.', ...
            nameStr, numel(v), Nfaces);
        v = v(1:Nfaces);
    elseif numel(v) < Nfaces
        warning('%s has %d entries, expected %d. Padding zeros.', ...
            nameStr, numel(v), Nfaces);
        v = [v; zeros(Nfaces - numel(v),1)];
    end
end

function plot_linear_patch(X,Y,Z,dataVals,ttl)
    dataVals = dataVals(:);
    dataVals(~isfinite(dataVals)) = 0;

    if numel(dataVals) ~= size(X,1)
        error('plot_linear_patch: %d values for %d faces.', numel(dataVals), size(X,1));
    end

    patch(X.', Y.', Z.', abs(dataVals), ...
        'FaceAlpha',1, ...
        'EdgeAlpha',0.2);

    title(ttl,'Interpreter','none');
    xlabel('X [m]');
    ylabel('Y [m]');
    zlabel('Z [m]');
    colorbar('eastoutside');
    axis equal tight;
    view(30,30);
end

function plot_log_patch(X,Y,Z,dataVals,ttl)
    dataVals = dataVals(:);
    dataVals(~isfinite(dataVals)) = 0;

    if numel(dataVals) ~= size(X,1)
        error('plot_log_patch: %d values for %d faces.', numel(dataVals), size(X,1));
    end

    patch(X.', Y.', Z.', log10(abs(dataVals) + eps), ...
        'FaceAlpha',1, ...
        'EdgeAlpha',0.2);

    title(ttl,'Interpreter','none');
    xlabel('X [m]');
    ylabel('Y [m]');
    zlabel('Z [m]');
    colorbar('eastoutside');
    axis equal tight;
    view(30,30);
end