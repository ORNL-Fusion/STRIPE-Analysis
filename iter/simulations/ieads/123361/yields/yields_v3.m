% Clear workspace and close all figures
% clear; clc; close all;

%% ======================= GEOMETRY ==========================
if ~exist('x1', 'var')
    fid = fopen('gitrGeometryPointPlane3d.cfg');
    if fid == -1
        error('Failed to open gitrGeometryPointPlane3d.cfg');
    end

    for i = 1:20
        tline = fgetl(fid);
        if i > 2 && ischar(tline)
            evalc(tline);
        end
    end
    fclose(fid);
end

subset = 1:length(x1);
X = [transpose(x1(subset)), transpose(x2(subset)), transpose(x3(subset))];
Y = [transpose(y1(subset)), transpose(y2(subset)), transpose(y3(subset))];
Z = [transpose(z1(subset)), transpose(z2(subset)), transpose(z3(subset))];

if ~exist('X', 'var') || ~exist('Y', 'var') || ~exist('Z', 'var')
    error('Geometry variables X, Y, and Z are not defined. Check gitrGeometryPointPlane3d.cfg.');
end

Nfaces = size(X,1);

%% ======================= PLOTTING HELPERS =============================
function plot_linear_patch(X,Y,Z,dataVals,ttl)
    dataVals = dataVals(:);
    if numel(dataVals) ~= size(X,1)
        error('plot_linear_patch: C-size mismatch: %d values for %d faces', ...
              numel(dataVals), size(X,1));
    end
    patch(transpose(X), transpose(Y), transpose(Z), abs(dataVals), ...
          'FaceAlpha', 1, 'EdgeAlpha', 0.3);
    title(ttl,'Interpreter','none');
    xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]');
    colorbar('eastoutside');
    axis equal tight;
    view(30,30);
end

function plot_log_patch(X,Y,Z,dataVals,ttl)
    dataVals = dataVals(:);
    if numel(dataVals) ~= size(X,1)
        error('plot_log_patch: C-size mismatch: %d values for %d faces', ...
              numel(dataVals), size(X,1));
    end
    patch(transpose(X), transpose(Y), transpose(Z), log10(abs(dataVals)+eps), ...
          'FaceAlpha', 1, 'EdgeAlpha', 0.3);
    title(ttl,'Interpreter','none');
    xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]');
    colorbar('eastoutside');
    axis equal tight;
    view(30,30);
end

%% ======================= FILE LISTS ==============================
labels = {'D+','Ne1+','Ne2+','Ne3+','Ne4+','Ne5+','Ne6+','Ne7+','Ne8+','Ne9+','Ne10+'};

target_files = {
    '../ieads_D+/Targets_D+.txt', ...
    '../ieads_Ne1+/Targets_Ne1+.txt', ...
    '../ieads_Ne2+/Targets_Ne2+.txt', ...
    '../ieads_Ne3+/Targets_Ne3+.txt', ...
    '../ieads_Ne4+/Targets_Ne4+.txt', ...
    '../ieads_Ne5+/Targets_Ne5+.txt', ...
    '../ieads_Ne6+/Targets_Ne6+.txt', ...
    '../ieads_Ne7+/Targets_Ne7+.txt', ...
    '../ieads_Ne8+/Targets_Ne8+.txt', ...
    '../ieads_Ne9+/Targets_Ne9+.txt', ...
    '../ieads_Ne10+/Targets_Ne10+.txt'};

yield_files = {
    '../ieads_D+/yields_D+.csv', ...
    '../ieads_Ne1+/yields_Ne1+.csv', ...
    '../ieads_Ne2+/yields_Ne2+.csv', ...
    '../ieads_Ne3+/yields_Ne3+.csv', ...
    '../ieads_Ne4+/yields_Ne4+.csv', ...
    '../ieads_Ne5+/yields_Ne5+.csv', ...
    '../ieads_Ne6+/yields_Ne6+.csv', ...
    '../ieads_Ne7+/yields_Ne7+.csv', ...
    '../ieads_Ne8+/yields_Ne8+.csv', ...
    '../ieads_Ne9+/yields_Ne9+.csv', ...
    '../ieads_Ne10+/yields_Ne10+.csv'};

%% ======================= INITIALIZE TOTALS ============================
total_ero_flux   = zeros(Nfaces,1);   % sum_s (n_s * Y_s * v_s)
total_ion_flux   = zeros(Nfaces,1);   % sum_s (n_s * v_s)
density_total    = zeros(Nfaces,1);   % sum_s (n_s)
yield_num_total  = zeros(Nfaces,1);   % sum_s (Y_s * n_s * v_s)

% NEW: direct sum of per-species yields
yield_total_sum  = zeros(Nfaces,1);   % sum_s (Y_s)

% store for optional effective-yield plots
all_yield = zeros(Nfaces, numel(labels));
all_flux  = zeros(Nfaces, numel(labels));
all_nvf   = zeros(Nfaces, numel(labels));   % n*v per species

%% ======================= PER-SPECIES LOOP =============================
for i = 1:numel(labels)
    lbl   = labels{i};
    data  = readmatrix(target_files{i});
    Yfile = readmatrix(yield_files{i});

    % ---- choose correct charge column ----
    if i == 1
        chargeIdx = 1;   % D+
    else
        tok = regexp(lbl,'Ne(\d+)\+','tokens','once');
        if isempty(tok)
            error('Could not parse charge from %s', lbl);
        end
        chargeIdx = str2double(tok{1});
    end

    if chargeIdx < 1 || chargeIdx > size(Yfile,2)
        error('Yields file %s has only %d columns; need column %d for %s', ...
              yield_files{i}, size(Yfile,2), chargeIdx, lbl);
    end

    % ---- align yield array with Nfaces ----
    ycol = Yfile(:,chargeIdx);

    if numel(ycol) == Nfaces+1 && ycol(1) == 0
        Y_align = ycol;
    else
        Y_align = [0; ycol];
    end

    if numel(Y_align) > Nfaces
        Y_align = Y_align(1:Nfaces);
    elseif numel(Y_align) < Nfaces
        error('After padding, yields len=%d < faces=%d for %s.', ...
              numel(Y_align), Nfaces, lbl);
    end

    % ---- density and velocity ----
    v_data = data(:,5);

    if i == 1
        dens_data = data(:,2);   % D+ density
    else
        dens_data = data(:,11);  % Ne density
    end

    if numel(dens_data) < Nfaces || numel(v_data) < Nfaces
        error('Density/flow shorter than Nfaces for %s (dens=%d, flow=%d, faces=%d).', ...
              lbl, numel(dens_data), numel(v_data), Nfaces);
    end

    dens_data = dens_data(1:Nfaces);
    v_data    = v_data(1:Nfaces);

    % ---- ion flux and erosion flux ----
    ion_flux_species = dens_data(:) .* v_data(:);          % n_s * v_s
    ero_flux_species = ion_flux_species .* Y_align(:);     % n_s * v_s * Y_s

    % ---- accumulate totals directly ----
    total_ion_flux  = total_ion_flux  + ion_flux_species;
    total_ero_flux  = total_ero_flux  + ero_flux_species;
    density_total   = density_total   + dens_data(:);
    yield_num_total = yield_num_total + ero_flux_species;

    % NEW: direct total yield from yield data only
    yield_total_sum = yield_total_sum + Y_align(:);

    % ---- store optional arrays ----
    all_yield(:,i) = Y_align(:);
    all_flux(:,i)  = ero_flux_species;
    all_nvf(:,i)   = ion_flux_species;

    % ---- per-species plots ----
    figure('Name', ['Species: ', lbl]);

    subplot(2,1,1);
    plot_linear_patch(X,Y,Z,Y_align, ['Sputtering Yield for ', lbl]);
    set(gca,'FontSize',10);

    subplot(2,1,2);
    plot_log_patch(X,Y,Z,ero_flux_species, ['Gross Erosion Flux for ', lbl]);
    set(gca,'FontSize',10);
end

%% ======================= EFFECTIVE TOTAL YIELD ========================
% Flux-weighted effective yield:
% Y_eff,total = total_ero_flux / total_ion_flux
yield_eff_total = nan(Nfaces,1);
maskTot = total_ion_flux ~= 0;
yield_eff_total(maskTot) = total_ero_flux(maskTot) ./ total_ion_flux(maskTot);

%% ======================= D AND NE EFFECTIVE YIELDS ====================
% D+
yield_eff_D = all_yield(:,1);

% Neon combined effective yield
ne_idx = 2:size(all_yield,2);
yield_eff_Ne = nan(Nfaces,1);

num_Ne = sum(all_yield(:,ne_idx) .* all_nvf(:,ne_idx), 2, 'omitnan');
den_Ne = sum(all_nvf(:,ne_idx), 2, 'omitnan');

maskNe = den_Ne ~= 0;
yield_eff_Ne(maskNe) = num_Ne(maskNe) ./ den_Ne(maskNe);

%% ======================= TOTAL EROSION FLUX PLOTS =====================
figure('Name','Total Erosion Flux (All species)');

subplot(2,1,1);
plot_linear_patch(X,Y,Z,total_ero_flux, 'Total Erosion Flux (linear)');
set(gca,'FontSize',10);

subplot(2,1,2);
plot_log_patch(X,Y,Z,total_ero_flux, 'Total Erosion Flux (log10)');
set(gca,'FontSize',10);

%% ======================= TOTAL YIELD SUM PLOTS ========================
figure('Name','Total Yield from Yield Data Only');

subplot(2,1,1);
plot_linear_patch(X,Y,Z,yield_total_sum, 'Total Yield = Sum of Per-Species Yields');
set(gca,'FontSize',10);

subplot(2,1,2);
plot_log_patch(X,Y,Z,yield_total_sum, 'Total Yield = Sum of Per-Species Yields (log10)');
set(gca,'FontSize',10);

%% ======================= EFFECTIVE YIELD PLOTS ========================
figure('Name','Effective Yields','Units','normalized','Position',[0.1 0.1 0.9 0.45]);

subplot(1,4,1);
patch(transpose(X), transpose(Y), transpose(Z), yield_eff_D, ...
      'FaceAlpha',1,'EdgeAlpha',0.3);
title('Effective Yield — D^+','Interpreter','none');
xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]');
colorbar('eastoutside');
axis equal tight;
view(30,30);
set(gca,'FontSize',10);

subplot(1,4,2);
patch(transpose(X), transpose(Y), transpose(Z), yield_eff_Ne, ...
      'FaceAlpha',1,'EdgeAlpha',0.3);
title('Effective Yield — Ne (all charge states)','Interpreter','none');
xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]');
colorbar('eastoutside');
axis equal tight;
view(30,30);
set(gca,'FontSize',10);

subplot(1,4,3);
patch(transpose(X), transpose(Y), transpose(Z), yield_eff_total, ...
      'FaceAlpha',1,'EdgeAlpha',0.3);
title('Effective Yield — Total Flux Weighted','Interpreter','none');
xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]');
colorbar('eastoutside');
axis equal tight;
view(30,30);
set(gca,'FontSize',10);

subplot(1,4,4);
patch(transpose(X), transpose(Y), transpose(Z), yield_total_sum, ...
      'FaceAlpha',1,'EdgeAlpha',0.3);
title('Total Yield — Direct Sum of Species Yields','Interpreter','none');
xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]');
colorbar('eastoutside');
axis equal tight;
view(30,30);
set(gca,'FontSize',10);

%% ======================= OPTIONAL MATCHED COLOR LIMITS ================
allEff = [yield_eff_D(:); yield_eff_Ne(:); yield_eff_total(:); yield_total_sum(:)];
finiteVals = allEff(isfinite(allEff));
if ~isempty(finiteVals)
    cmin = min(finiteVals);
    cmax = max(finiteVals);
    subplot(1,4,1); caxis([cmin cmax]);
    subplot(1,4,2); caxis([cmin cmax]);
    subplot(1,4,3); caxis([cmin cmax]);
    subplot(1,4,4); caxis([cmin cmax]);
end

%% ======================= WRITE OUTPUT FILES ===========================
writematrix(total_ero_flux,   'total_erosion_flux.txt');
writematrix(total_ion_flux,   'total_ion_flux.txt');
writematrix(yield_eff_total,  'total_effective_yield.txt');
writematrix(yield_total_sum,  'total_yield_sum_from_species.txt');
writematrix(density_total,    'density_total_all_species.txt');