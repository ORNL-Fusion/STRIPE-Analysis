% total_yields_iter.m
% Clear workspace and close all figures
% clear; clc; close all;

%% ======================= GEOMETRY ==========================
% Load geometry if not already in workspace
if ~exist('x1', 'var')
    fid = fopen('gitrGeometryPointPlane3d.cfg');
    if fid == -1
        error('Failed to open gitrGeometryPointPlane3d.cfg');
    end
    % File has two header lines, then lines that define x1..z3 etc.
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

% Yields
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

%% ======================= YIELD LOADER (ROBUST) =========================
function y = load_yield_for_charge(yieldFile, Nfaces, chargeCol)
    % Robust loader that tolerates off-by-one (or small) mismatches by pad/trim.
    % Supports:
    %  (1) nr x 1 or nr x K (face-aligned)
    %  (2) M x (K+1) with first col = face index (1..Nfaces)
    % If face-aligned length differs from Nfaces, it pads with zeros or trims.

    M = readmatrix(yieldFile);

    if isempty(M)
        error('Yield file is empty: %s', yieldFile);
    end

    [nr,nc] = size(M);

    % --- Case B: first column is an index column ---
    if nc >= 2
        idx = M(:,1);
        if all(isfinite(idx)) && all(idx == round(idx)) && all(idx >= 1) && all(idx <= Nfaces)
            y = zeros(Nfaces,1);
            c = min(max(2, chargeCol+1), nc); % +1 because col1 is index
            vals = M(:,c);
            vals(~isfinite(vals)) = 0;
            y(idx) = vals;
            return;
        end
    end

    % --- Case A: face-aligned (no explicit indices) ---
    if nc == 1
        yraw = M(:,1);
    else
        c = min(max(1, chargeCol), nc);
        yraw = M(:,c);
    end

    % Clean non-finite
    yraw(~isfinite(yraw)) = 0;

    % If exact match, done
    if numel(yraw) == Nfaces
        y = yraw(:);
        return;
    end

    % Tolerate mismatch: pad or trim
    y = zeros(Nfaces,1);
    ncopy = min(Nfaces, numel(yraw));
    y(1:ncopy) = yraw(1:ncopy);

    if numel(yraw) < Nfaces
        warning('Yield length mismatch for %s: got %d, expected %d. Padded with zeros at end.', ...
            yieldFile, numel(yraw), Nfaces);
    else
        warning('Yield length mismatch for %s: got %d, expected %d. Trimmed extra rows.', ...
            yieldFile, numel(yraw), Nfaces);
    end
end

%% ======================= FIGURE: D+ vs TOTAL NEON YIELDS ===============
% D+ yield (chargeCol=1)
yD = load_yield_for_charge(yield_files{1}, Nfaces, 1);

% Total neon yield = sum over Ne1..Ne10 using matching charge column
yNeTot = zeros(Nfaces,1);
for q = 1:10
    % yield_files{q+1} corresponds to Ne(q)+
    yq = load_yield_for_charge(yield_files{q+1}, Nfaces, q);
    yNeTot = yNeTot + yq;
end

figure('Name','Yields: D+ and Total Neon','Color','w');
subplot(1,2,1);
plot_linear_patch(X,Y,Z,yD,'D+ yields');
view(-30,30);      % exact opposite of view(30,30)

subplot(1,2,2);
plot_linear_patch(X,Y,Z,yNeTot,'Total Neon yields (Ne1+..Ne10+)');

% If you prefer log scale, swap plot_linear_patch -> plot_log_patch above.
view(-30,30);      % exact opposite of view(30,30)    % mirror left–right (most common fix)