clear; clc; close all;

%% ======================= GEOMETRY ==========================
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

subset = 1:length(x1);

X = [x1(subset)', x2(subset)', x3(subset)'];
Y = [y1(subset)', y2(subset)', y3(subset)'];
Z = [z1(subset)', z2(subset)', z3(subset)'];

Nfaces = size(X,1);

%% ======================= FACE AREA =========================
A = 0.5 * sqrt( ...
    ((Y(:,2)-Y(:,1)).*(Z(:,3)-Z(:,1)) - (Z(:,2)-Z(:,1)).*(Y(:,3)-Y(:,1))).^2 + ...
    ((Z(:,2)-Z(:,1)).*(X(:,3)-X(:,1)) - (X(:,2)-X(:,1)).*(Z(:,3)-Z(:,1))).^2 + ...
    ((X(:,2)-X(:,1)).*(Y(:,3)-Y(:,1)) - (Y(:,2)-Y(:,1)).*(X(:,3)-X(:,1))).^2 );

%% ======================= PLOT HELPERS ======================
function plot_linear_patch(X,Y,Z,dataVals,ttl)
    patch(X',Y',Z',dataVals(:),'FaceAlpha',1,'EdgeAlpha',0.3);
    title(ttl,'Interpreter','none');
    xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]');
    colorbar; axis equal tight; view(-30,30);
end

function plot_log_patch(X,Y,Z,dataVals,ttl)
    patch(X',Y',Z',log10(abs(dataVals(:))+eps),'FaceAlpha',1,'EdgeAlpha',0.3);
    title(ttl,'Interpreter','none');
    xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]');
    colorbar; axis equal tight; view(-30,30);
end

%% ======================= FILE LISTS ========================
labels = {'D+','Ne1+','Ne2+','Ne3+','Ne4+','Ne5+', ...
          'Ne6+','Ne7+','Ne8+','Ne9+','Ne10+'};

target_files = {
    '../ieads_D+/Targets_D+.txt'
    '../ieads_Ne1+/Targets_Ne1+.txt'
    '../ieads_Ne2+/Targets_Ne2+.txt'
    '../ieads_Ne3+/Targets_Ne3+.txt'
    '../ieads_Ne4+/Targets_Ne4+.txt'
    '../ieads_Ne5+/Targets_Ne5+.txt'
    '../ieads_Ne6+/Targets_Ne6+.txt'
    '../ieads_Ne7+/Targets_Ne7+.txt'
    '../ieads_Ne8+/Targets_Ne8+.txt'
    '../ieads_Ne9+/Targets_Ne9+.txt'
    '../ieads_Ne10+/Targets_Ne10+.txt'};

yield_files = {
    '../ieads_D+/yields_D+.csv'
    '../ieads_Ne1+/yields_Ne1+.csv'
    '../ieads_Ne2+/yields_Ne2+.csv'
    '../ieads_Ne3+/yields_Ne3+.csv'
    '../ieads_Ne4+/yields_Ne4+.csv'
    '../ieads_Ne5+/yields_Ne5+.csv'
    '../ieads_Ne6+/yields_Ne6+.csv'
    '../ieads_Ne7+/yields_Ne7+.csv'
    '../ieads_Ne8+/yields_Ne8+.csv'
    '../ieads_Ne9+/yields_Ne9+.csv'
    '../ieads_Ne10+/yields_Ne10+.csv'};

%% ======================= LOADERS ===========================
function Mout = load_face_matrix(fname,Nfaces)
    M = readmatrix(fname);
    if isempty(M)
        error('Empty file: %s',fname);
    end

    if size(M,1) >= Nfaces
        Mout = M(1:Nfaces,:);
    else
        Mout = [M; zeros(Nfaces-size(M,1),size(M,2))];
        warning('Padded %s with zeros.',fname);
    end

    Mout(~isfinite(Mout)) = 0;
end

function y = load_yield_for_charge(fname,Nfaces,chargeCol)
    M = readmatrix(fname);
    if isempty(M)
        error('Empty yield file: %s',fname);
    end

    [~,nc] = size(M);

    if nc >= 2
        idx = M(:,1);
        if all(isfinite(idx)) && all(idx == round(idx)) && ...
           all(idx >= 1) && all(idx <= Nfaces)

            y = zeros(Nfaces,1);
            c = min(chargeCol+1,nc);
            vals = M(:,c);
            vals(~isfinite(vals)) = 0;
            y(idx) = vals;
            return;
        end
    end

    c = min(chargeCol,nc);
    yraw = M(:,c);
    yraw(~isfinite(yraw)) = 0;

    y = zeros(Nfaces,1);
    ncopy = min(Nfaces,numel(yraw));
    y(1:ncopy) = yraw(1:ncopy);
end

%% ======================= LOAD TARGET DATA ==================
data = cell(numel(labels),1);
for i = 1:numel(labels)
    data{i} = load_face_matrix(target_files{i},Nfaces);
end

%% ======================= DENSITY AND FLOW ==================
nList = cell(numel(labels),1);
vList = cell(numel(labels),1);

nList{1} = data{1}(:,2);
vList{1} = data{1}(:,5);

for i = 2:numel(labels)
    nList{i} = data{i}(:,11);
    vList{i} = data{i}(:,5);
end

%% ======================= VELOCITY DIAGNOSTIC ===============
disp(' ');
disp('Velocity diagnostic:');
for i = 1:numel(labels)
    fprintf('%6s : min(v)=%.3e, max(v)=%.3e, positive=%d, negative=%d\n', ...
        labels{i}, min(vList{i}), max(vList{i}), ...
        sum(vList{i}>0), sum(vList{i}<0));
end

%% ======================= YIELD MAPS ========================
yD = load_yield_for_charge(yield_files{1},Nfaces,1);

yNeTot = zeros(Nfaces,1);
for q = 1:10
    yNeTot = yNeTot + load_yield_for_charge(yield_files{q+1},Nfaces,q);
end

figure('Name','Yields: D+ and Total Neon','Color','w');

subplot(1,2,1);
plot_linear_patch(X,Y,Z,yD,'D+ yields');

subplot(1,2,2);
plot_linear_patch(X,Y,Z,yNeTot,'Total Neon yields');

%% ======================= EROSION CALCULATION ===============
erosion_rate_per_ion = zeros(numel(labels),1);
facewise_erosion_flux_density = zeros(Nfaces,numel(labels));
facewise_erosion_rate = zeros(Nfaces,numel(labels));

for i = 1:numel(labels)

    if i == 1
        chargeCol = 1;
    else
        chargeCol = i - 1;
    end

    n_i = nList{i}(:);
    v_i = vList{i}(:);
    y_i = load_yield_for_charge(yield_files{i},Nfaces,chargeCol);

    % Positive velocity is treated as incoming.
    Gamma_in = n_i .* max(v_i,0);

    erosion_flux_density = Gamma_in .* y_i;
    erosion_rate_face = erosion_flux_density .* A;

    facewise_erosion_flux_density(:,i) = erosion_flux_density;
    facewise_erosion_rate(:,i) = erosion_rate_face;

    erosion_rate_per_ion(i) = sum(erosion_rate_face,'omitnan');
end

%% ======================= BAR CHART =========================
figure('Name','Area-Integrated Gross Erosion Rate by Ion','Color','w');
bar(erosion_rate_per_ion);
set(gca,'XTick',1:numel(labels), ...
        'XTickLabel',labels, ...
        'XTickLabelRotation',45);
xlabel('Ion species');
ylabel('\Sigma \Gamma_{in} Y A');
title('Area-integrated gross erosion contribution');
grid on;

%% ======================= LOG BAR CHART =====================
figure('Name','Log Area-Integrated Gross Erosion Rate by Ion','Color','w');
bar(erosion_rate_per_ion);
set(gca,'YScale','log');
set(gca,'XTick',1:numel(labels), ...
        'XTickLabel',labels, ...
        'XTickLabelRotation',45);
xlabel('Ion species');
ylabel('\Sigma \Gamma_{in} Y A');
title('Area-integrated gross erosion contribution');
grid on;

%% ======================= TOTAL FLUX MAP ====================
total_flux_density = sum(facewise_erosion_flux_density,2,'omitnan');

figure('Name','Total Gross Erosion Flux Density Map','Color','w');
plot_linear_patch(X,Y,Z,total_flux_density, ...
    'Total gross erosion flux density');

figure('Name','Log Total Gross Erosion Flux Density Map','Color','w');
plot_log_patch(X,Y,Z,total_flux_density, ...
    'Log total gross erosion flux density');

%% ======================= TOTAL RATE MAP ====================
total_rate_face = sum(facewise_erosion_rate,2,'omitnan');

figure('Name','Area-Integrated Erosion Rate per Face','Color','w');
plot_linear_patch(X,Y,Z,total_rate_face, ...
    'Area-integrated erosion rate per face');

%% ======================= PRINT VALUES ======================
disp(' ');
disp('Area-integrated gross erosion rate by ion species:');

for i = 1:numel(labels)
    fprintf('%6s : %.6e\n',labels{i},erosion_rate_per_ion(i));
end

fprintf('TOTAL  : %.6e\n',sum(erosion_rate_per_ion,'omitnan'));