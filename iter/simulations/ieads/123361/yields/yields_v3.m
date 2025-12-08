% Clear workspace and close all figures
% clear; clc; close all;

%% ======================= GEOMETRY ==========================
% Load geometry if not already in workspace
if ~exist('x1', 'var')
    fid = fopen('gitrGeometryPointPlane3d.cfg');
    if fid == -1
        error('Failed to open gitrGeometryPointPlane3d.cfg');
    end
    % File has two header lines, then 18 lines that define x1..z3 etc.
    for i = 1:20
        tline = fgetl(fid);
        if i > 2 && ischar(tline)
            evalc(tline); %#ok<EVLC>
        end
    end
    fclose(fid);
end

% Define subset of faces and build X,Y,Z from geometry
subset = 1:length(x1);              % Nfaces
% Store as Nfaces x 3, we will transpose when calling patch (to 3 x Nfaces)
X = [transpose(x1(subset)), transpose(x2(subset)), transpose(x3(subset))];
Y = [transpose(y1(subset)), transpose(y2(subset)), transpose(y3(subset))];
Z = [transpose(z1(subset)), transpose(z2(subset)), transpose(z3(subset))];

% Ensure X, Y, Z are defined
if ~exist('X', 'var') || ~exist('Y', 'var') || ~exist('Z', 'var')
    error('Geometry variables X, Y, and Z are not defined. Check gitrGeometryPointPlane3d.cfg.');
end

Nfaces = size(X,1);

%% ======================= PLOTTING HELPERS =============================
function plot_linear_patch(X,Y,Z,dataVals,ttl)
    % Expects X/Y/Z as Nfacesx3; transpose to 3xNfaces for patch
    dataVals = dataVals(:);
    if numel(dataVals) ~= size(X,1)
        error('plot_linear_patch: C-size mismatch: %d values for %d faces', numel(dataVals), size(X,1));
    end
    patch(transpose(X), transpose(Y), transpose(Z), abs(dataVals), ...
          'FaceAlpha', 1, 'EdgeAlpha', 0.3);
    title(ttl,'Interpreter','none');
    xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]');
    colorbar('eastoutside'); axis equal tight; view(30,30);
end

function plot_log_patch(X,Y,Z,dataVals,ttl)
    dataVals = dataVals(:);
    if numel(dataVals) ~= size(X,1)
        error('plot_log_patch: C-size mismatch: %d values for %d faces', numel(dataVals), size(X,1));
    end
    patch(transpose(X), transpose(Y), transpose(Z), log10(abs(dataVals)+eps), ...
          'FaceAlpha', 1, 'EdgeAlpha', 0.3);
    title(ttl,'Interpreter','none');
    xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]');
    colorbar('eastoutside'); axis equal tight; view(30,30);
end

%% ======================= FILE LISTS (D+ + Ne1..Ne10) ==================
labels = {'D+','Ne1+','Ne2+','Ne3+','Ne4+','Ne5+','Ne6+','Ne7+','Ne8+','Ne9+','Ne10+'};

% Targets
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

% Yields (multi-column allowed; we will pick the correct charge column)
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

%% ======================= PER-SPECIES CALC & PLOTS =====================
all_flux  = zeros(Nfaces, numel(labels));
all_dens  = zeros(Nfaces, numel(labels));
all_flow  = zeros(Nfaces, numel(labels));
all_yield = zeros(Nfaces, numel(labels));

for i = 1:numel(labels)
    lbl = labels{i};
    data = readmatrix(target_files{i});
    Yfile = readmatrix(yield_files{i});

    % ---- choose charge column ----
    if i == 1
        chargeIdx = 1;   % D+
    else
        tok = regexp(lbl,'Ne(\d+)\+','tokens','once');
        if isempty(tok), error('Could not parse charge from %s', lbl); end
        chargeIdx = str2double(tok{1});
    end
    if chargeIdx < 1 || chargeIdx > size(Yfile,2)
        error('Yields file %s has only %d columns; need column %d for %s', ...
              yield_files{i}, size(Yfile,2), chargeIdx, lbl);
    end

    % ---- convention: ALWAYS use [0; yields(:,charge)] ----
    ycol = Yfile(:,chargeIdx);
    if numel(ycol) == Nfaces+1 && ycol(1) == 0
        Y_align = ycol;                 % already padded
    else
        Y_align = [0; ycol];            % prepend leading 0
    end
    % Make sure length matches Nfaces
    if numel(Y_align) > Nfaces
        Y_align = Y_align(1:Nfaces);
    elseif numel(Y_align) < Nfaces
        error('After padding, yields len=%d < faces=%d for %s.', numel(Y_align), Nfaces, lbl);
    end

    % ---- Density & flow selection ----
    v_data  = data(:,5);                 % flow = col 5
    if i == 1
        dens_data = data(:,2);           % D+ density = col 2
    else
        dens_data = data(:,11);          % Ne density = col 11
    end
    if numel(dens_data) < Nfaces || numel(v_data) < Nfaces
        error('Density/flow shorter than Nfaces for %s (dens=%d, flow=%d, faces=%d).', ...
              lbl, numel(dens_data), numel(v_data), Nfaces);
    end
    dens_data = dens_data(1:Nfaces);
    v_data    = v_data(1:Nfaces);

    % ---- Erosion flux (per face) ----
    ero_data = Y_align(:) .* dens_data(:) .* v_data(:);
    all_flux(:,i)  = ero_data(:);
all_dens(:,i)  = dens_data(:);
all_flow(:,i)  = v_data(:);
all_yield(:,i) = Y_align(:);

    % ---- Plots per species ----
    figure('Name', ['Species: ', lbl]);
    subplot(2,1,1);
    plot_linear_patch(X,Y,Z,Y_align, ['Sputtering Yield for ', lbl]);
    set(gca, 'FontSize', 10);

    subplot(2,1,2);
    plot_log_patch(X,Y,Z,ero_data, ['Gross Erosion Flux for ', lbl]);
    set(gca, 'FontSize', 10);
    clim([14 20])
end

%% ======================= TOTAL (SUM OVER SPECIES) =====================
flux_total    = sum(all_flux,  2, 'omitnan');
density_total = sum(all_dens,  2, 'omitnan');
yields_total  = sum(all_yield, 2, 'omitnan');
flow_total    = mean(all_flow, 2, 'omitnan');  % simple average across species

figure('Name','Total Erosion Flux (All species)');
subplot(2,1,1);
plot_linear_patch(X,Y,Z,flux_total, 'Total Erosion Flux (linear)');
set(gca, 'FontSize', 10);

subplot(2,1,2);
plot_log_patch(X,Y,Z,flux_total, 'Total Erosion Flux (log10)');
set(gca, 'FontSize', 10);

% ============== WRITE TOTALS TO DISK ========================
writematrix(flux_total,     'flux_total_all_species.txt');
writematrix(density_total,  'density_total_all_species.txt');
writematrix(yields_total,   'yields_total_all_species.txt');
writematrix(flow_total,     'flow_total_all_species.txt');
