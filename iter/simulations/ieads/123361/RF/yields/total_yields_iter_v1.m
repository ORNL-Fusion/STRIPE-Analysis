% total_yields_iter.m
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

%% ======================= FACE AREA =========================
% Triangle area from 3 vertices:
% area = 0.5 * | (p2-p1) x (p3-p1) |

P1 = [X(:,1), Y(:,1), Z(:,1)];
P2 = [X(:,2), Y(:,2), Z(:,2)];
P3 = [X(:,3), Y(:,3), Z(:,3)];

v1 = P2 - P1;
v2 = P3 - P1;
cp = cross(v1, v2, 2);

area = 0.5 * sqrt(sum(cp.^2, 2));   % [m^2], one value per face
area(~isfinite(area)) = 0;

%% ======================= PLOTTING HELPERS =============================
function plot_linear_patch(X,Y,Z,dataVals,ttl)
    dataVals = dataVals(:);
    if numel(dataVals) ~= size(X,1)
        error('plot_linear_patch: C-size mismatch: %d values for %d faces', ...
            numel(dataVals), size(X,1));
    end

    patch(transpose(X), transpose(Y), transpose(Z), abs(dataVals), ...
        'FaceColor', 'flat', ...
        'FaceAlpha', 1, ...
        'EdgeAlpha', 0.15);

    title(ttl,'Interpreter','none');
    xlabel('X [m]');
    ylabel('Y [m]');
    zlabel('Z [m]');
    colorbar('eastoutside');
    axis equal tight;
    view(-30,30);
    % lighting gouraud;
    % camlight headlight;
end

function plot_log_patch(X,Y,Z,dataVals,ttl)
    dataVals = dataVals(:);
    if numel(dataVals) ~= size(X,1)
        error('plot_log_patch: C-size mismatch: %d values for %d faces', ...
            numel(dataVals), size(X,1));
    end

    patch(transpose(X), transpose(Y), transpose(Z), log10(abs(dataVals)+eps), ...
        'FaceColor', 'flat', ...
        'FaceAlpha', 1, ...
        'EdgeAlpha', 0.15);

    title(ttl,'Interpreter','none');
    xlabel('X [m]');
    ylabel('Y [m]');
    zlabel('Z [m]');
    colorbar('eastoutside');
    axis equal tight;
    view(-30,30);
    lighting gouraud;
    camlight headlight;
end

%% ======================= FILE LISTS =======================
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

%% ======================= ROBUST MATRIX LOADER ==========================
function Mout = load_face_matrix(fname, Nfaces)
    M = readmatrix(fname);

    if isempty(M)
        error('File is empty: %s', fname);
    end

    [nr,~] = size(M);

    if nr >= Nfaces
        Mout = M(1:Nfaces,:);
    else
        nMissing = Nfaces - nr;
        Mout = [M; zeros(nMissing, size(M,2))];
        warning('File %s has %d rows, expected %d. Padded with zeros.', ...
            fname, nr, Nfaces);
    end

    Mout(~isfinite(Mout)) = 0;
end

%% ======================= YIELD LOADER ======================
function y = load_yield_for_charge(yieldFile, Nfaces, chargeCol)
    M = readmatrix(yieldFile);

    if isempty(M)
        error('Yield file is empty: %s', yieldFile);
    end

    [~,nc] = size(M);

    % Case 1: indexed file
    if nc >= 2
        idx = M(:,1);
        if all(isfinite(idx)) && all(idx == round(idx)) && ...
                all(idx >= 1) && all(idx <= Nfaces)

            y = zeros(Nfaces,1);
            c = min(max(2, chargeCol+1), nc);
            vals = M(:,c);
            vals(~isfinite(vals)) = 0;
            y(idx) = vals;
            return;
        end
    end

    % Case 2: face-aligned file
    if nc == 1
        yraw = M(:,1);
    else
        c = min(max(1, chargeCol), nc);
        yraw = M(:,c);
    end

    yraw(~isfinite(yraw)) = 0;

    y = zeros(Nfaces,1);
    ncopy = min(Nfaces, numel(yraw));
    y(1:ncopy) = yraw(1:ncopy);

    if numel(yraw) ~= Nfaces
        warning('Yield length mismatch for %s: got %d, expected %d.', ...
            yieldFile, numel(yraw), Nfaces);
    end
end

%% ======================= LOAD TARGET DATA ==================
data_D    = load_face_matrix(target_files{1},  Nfaces);
data_Ne1  = load_face_matrix(target_files{2},  Nfaces);
data_Ne2  = load_face_matrix(target_files{3},  Nfaces);
data_Ne3  = load_face_matrix(target_files{4},  Nfaces);
data_Ne4  = load_face_matrix(target_files{5},  Nfaces);
data_Ne5  = load_face_matrix(target_files{6},  Nfaces);
data_Ne6  = load_face_matrix(target_files{7},  Nfaces);
data_Ne7  = load_face_matrix(target_files{8},  Nfaces);
data_Ne8  = load_face_matrix(target_files{9},  Nfaces);
data_Ne9  = load_face_matrix(target_files{10}, Nfaces);
data_Ne10 = load_face_matrix(target_files{11}, Nfaces);

%% ======================= DENSITY AND FLOW ==================
% Adjust these columns if your file layout differs
n_D    = data_D(1:Nfaces,2);     v_D    = data_D(1:Nfaces,5);
n_Ne1  = data_Ne1(1:Nfaces,11);  v_Ne1  = data_Ne1(1:Nfaces,5);
n_Ne2  = data_Ne2(1:Nfaces,11);  v_Ne2  = data_Ne2(1:Nfaces,5);
n_Ne3  = data_Ne3(1:Nfaces,11);  v_Ne3  = data_Ne3(1:Nfaces,5);
n_Ne4  = data_Ne4(1:Nfaces,11);  v_Ne4  = data_Ne4(1:Nfaces,5);
n_Ne5  = data_Ne5(1:Nfaces,11);  v_Ne5  = data_Ne5(1:Nfaces,5);
n_Ne6  = data_Ne6(1:Nfaces,11);  v_Ne6  = data_Ne6(1:Nfaces,5);
n_Ne7  = data_Ne7(1:Nfaces,11);  v_Ne7  = data_Ne7(1:Nfaces,5);
n_Ne8  = data_Ne8(1:Nfaces,11);  v_Ne8  = data_Ne8(1:Nfaces,5);
n_Ne9  = data_Ne9(1:Nfaces,11);  v_Ne9  = data_Ne9(1:Nfaces,5);
n_Ne10 = data_Ne10(1:Nfaces,11); v_Ne10 = data_Ne10(1:Nfaces,5);

nList = {n_D,n_Ne1,n_Ne2,n_Ne3,n_Ne4,n_Ne5,n_Ne6,n_Ne7,n_Ne8,n_Ne9,n_Ne10};
vList = {v_D,v_Ne1,v_Ne2,v_Ne3,v_Ne4,v_Ne5,v_Ne6,v_Ne7,v_Ne8,v_Ne9,v_Ne10};

%% ======================= FIGURE: D+ vs TOTAL NEON YIELDS ===============
yD = load_yield_for_charge(yield_files{1}, Nfaces, 1);

yNeTot = zeros(Nfaces,1);
for q = 1:10
    yq = load_yield_for_charge(yield_files{q+1}, Nfaces, q);
    yNeTot = yNeTot + yq;
end

figure('Name','Yields: D+ and Total Neon','Color','w');

subplot(2,1,1);
plot_linear_patch(X,Y,Z,yD,'D+ yields');

subplot(2,1,2);
plot_linear_patch(X,Y,Z,yNeTot,'Total Neon yields (Ne1+..Ne10+)');

%% ======================= AREA-INTEGRATED EROSION =======================
% erosion_flux_i(face) = |n_i * v_i * Yeff_i|
% erosion_i(face)      = erosion_flux_i(face) * area(face)
%
% This gives an integrated erosion rate contribution per face.

abs_erosion_rate_per_ion = zeros(numel(labels),1);
facewise_abs_erosion_rate = zeros(Nfaces, numel(labels));

for i = 1:numel(labels)

    if i == 1
        chargeCol = 1;
    else
        chargeCol = i - 1;
    end

    n_i = nList{i}(:);
    v_i = vList{i}(:);
    y_eff_i = load_yield_for_charge(yield_files{i}, Nfaces, chargeCol);

    erosion_flux_i = abs(n_i .* v_i .* y_eff_i);   % [something / m^2 / s]
    erosion_i = erosion_flux_i .* area;            % [something / s]

    erosion_i(~isfinite(erosion_i)) = 0;

    facewise_abs_erosion_rate(:,i) = erosion_i;
    abs_erosion_rate_per_ion(i) = sum(erosion_i, 'omitnan');
end

%% ======================= BAR CHART =========================
figure('Name','Area-Integrated Erosion by Ion Species','Color','w');

totalErosion = sum(abs_erosion_rate_per_ion, 'omitnan');
pct = 100 * abs_erosion_rate_per_ion / max(totalErosion, eps);

bar(abs_erosion_rate_per_ion);
hold on;

set(gca, 'XTick', 1:numel(labels), ...
         'XTickLabel', labels, ...
         'XTickLabelRotation', 45, ...
         'FontName','Arial', ...
         'FontSize',12);

xlabel('Ion species');
ylabel('Area-integrated erosion rate');
title('Area-integrated erosion contribution from each ion species');
grid on;

for i = 1:numel(labels)
    text(i, abs_erosion_rate_per_ion(i), ...
        sprintf('%.1f%%', pct(i)), ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','bottom', ...
        'FontSize',10, ...
        'FontWeight','bold');
end

%% ======================= LOG BAR CHART =====================
figure('Name','Log Area-Integrated Erosion by Ion Species','Color','w');
bar(abs_erosion_rate_per_ion);
set(gca, 'YScale', 'log');
set(gca, 'XTick', 1:numel(labels), ...
         'XTickLabel', labels, ...
         'XTickLabelRotation', 45, ...
         'FontName','Arial', ...
         'FontSize',12);
xlabel('Ion species');
ylabel('Area-integrated erosion rate (log scale)');
title('Area-integrated erosion contribution from each ion species');
grid on;

%% ======================= TOTAL EROSION MAP =================
total_abs_erosion_all_ions = sum(facewise_abs_erosion_rate, 2, 'omitnan');

figure('Name','Total Area-Integrated Erosion Map (All Ions)','Color','w');
plot_linear_patch(X,Y,Z,total_abs_erosion_all_ions, ...
    'Total area-integrated erosion (all ions)');

%% ======================= LOG TOTAL EROSION MAP =============
figure('Name','Log Total Area-Integrated Erosion Map (All Ions)','Color','w');
plot_log_patch(X,Y,Z,total_abs_erosion_all_ions, ...
    'Log total area-integrated erosion (all ions)');

%% ======================= OPTIONAL INDIVIDUAL MAPS ==========
% for i = 1:numel(labels)
%     figure('Name',['Erosion Map - ' labels{i}], 'Color','w');
%     plot_linear_patch(X,Y,Z,facewise_abs_erosion_rate(:,i), ...
%         ['Area-integrated erosion - ' labels{i}]);
% end

%% ======================= OPTIONAL SURFACE SUBSET / CDF =====
% If you have a valid surface-face index list called surf_inds, this gives:
%   erosion_sub
%   erosion_sub_cdf
%   erosion_rate
%
% Example based on total erosion from all ions:
if exist('surf_inds', 'var') && ~isempty(surf_inds)
    erosion = total_abs_erosion_all_ions(:);

    erosion = erosion(surf_inds);
    erosion_inds = find(erosion);
    erosion_sub = erosion(erosion_inds);
    erosion_sub_cdf = cumsum(erosion_sub);

    if ~isempty(erosion_sub_cdf) && erosion_sub_cdf(end) > 0
        erosion_rate = erosion_sub_cdf(end);
        erosion_sub_cdf = erosion_sub_cdf ./ erosion_sub_cdf(end);
    else
        erosion_rate = 0;
        erosion_sub_cdf = [];
    end

    disp(' ');
    fprintf('Subset erosion rate = %.6e\n', erosion_rate);
end

%% ======================= PRINT VALUES ======================
disp(' ');
disp('Area-integrated erosion by ion species:');
for i = 1:numel(labels)
    fprintf('%6s : %.6e  (%.2f%%)\n', labels{i}, abs_erosion_rate_per_ion(i), pct(i));
end
fprintf('TOTAL  : %.6e\n', totalErosion);