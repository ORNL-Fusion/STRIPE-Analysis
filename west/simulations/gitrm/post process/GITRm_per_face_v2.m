clc; clear; close all; format long

%% --- Load GITRm per-face data ---
path1 = '../WEST_runs_test/gitrm-surface2.nc';
ncid1 = netcdf.open(path1,'NOWRITE');

coords      = netcdf.getVar(ncid1, netcdf.inqVarID(ncid1,'Coordinates'));
grossPlasma = squeeze(netcdf.getVar(ncid1, netcdf.inqVarID(ncid1,'PlasmaErosion')));
grossDepo   = squeeze(netcdf.getVar(ncid1, netcdf.inqVarID(ncid1,'grossDeposition')));
grossEroRaw = squeeze(netcdf.getVar(ncid1, netcdf.inqVarID(ncid1,'grossErosion')));
aveSpyl     = squeeze(netcdf.getVar(ncid1, 6));  % 37476x1 (average sputtering yield)
spylCounts  = squeeze(netcdf.getVar(ncid1, 7));  % 37476x1 (counts used for aveSpyl)

netcdf.close(ncid1);

% Number of triangles (from data vectors)
sz = size(grossEroRaw,1);

%% --- Unpack triangle vertices (fast, no loop) ---
% coords is 9*sz x 1 (or 1 x 9*sz). Reshape to [9 x sz] then transpose to [sz x 9].
coords_mat = reshape(coords, [9, sz]).';
x1 = coords_mat(:,1); y1 = coords_mat(:,2); z1 = coords_mat(:,3);
x2 = coords_mat(:,4); y2 = coords_mat(:,5); z2 = coords_mat(:,6);
x3 = coords_mat(:,7); y3 = coords_mat(:,8); z3 = coords_mat(:,9);

%% --- Triangle area (3D) ---
P1 = [x1 y1 z1];  P2 = [x2 y2 z2];  P3 = [x3 y3 z3];
v1 = P2 - P1;     v2 = P3 - P1;
cp = cross(v1, v2, 2);                 % row-wise cross
ElementArea = 0.5 * sqrt(sum(cp.^2,2));% sz x 1
% Guard against degenerate zero-area elements
ElementArea(ElementArea==0) = eps;

%% --- Erosion pre-factors and per-area normalization ---
% erosion_rate = 4.090184441246757e+17;  % #196154 (units depend on your tally definition)
erosion_rate= 1.062786872720573e+19;  % # WEST-RF
nP           = 5E3;                    % primary particles
erosionPP    = erosion_rate / nP;

% Normalize per face by area (-> per-area rates/fluxes)
grossPlasmaA   = erosionPP .* grossPlasma   ./ ElementArea;
grossDepoA     = abs(erosionPP .* grossDepo)     ./ ElementArea;
grossEroRawA   = erosionPP .* grossEroRaw   ./ ElementArea;

% If "gross erosion" should include plasma-driven component:
grossEroA = grossPlasmaA + grossEroRawA;   % per-area gross erosion
% grossEroA = grossEroRawA;   % per-area gross erosion
% grossEroA = grossPlasmaA;   % per-area gross erosion
netEroA   = grossEroA   - grossDepoA;      % per-area net erosion
netDepoA   = grossDepoA-grossEroRawA;      % per-area net erosion

%% =========================
%  Figure 1: 3D mesh subplots
%  ==========================
vars3D = { ...
    grossEroA,  '3D: gross erosion', 'particles/m^2/s'; ...
    grossDepoA, '3D: gross deposition',    'particles/m^2/s'; ...
    netEroA,    '3D: net erosion',   'particles/m^2/s' ...
    };

figure('Name','3D mesh (log-scale, per area)','Position',[100 100 1200 450]);

for k = 1:3
    subplot(1,3,k)
    C3D = (abs(vars3D{k,1}));
    C3D(~isfinite(C3D)) = NaN;

    patch( ...
        [x1 x2 x3]', ...
        [y1 y2 y3]', ...
        [z1 z2 z3]', ...
        C3D', ...
        'FaceColor','flat', ...
        'EdgeColor','none', ...
        'FaceAlpha',1, ...
        'EdgeAlpha',0.1);
    set(gca, 'ColorScale', 'log', 'FontSize', 16);

    axis equal
    xlabel('X'); ylabel('Y'); zlabel('Z');
    grid on; view(3);
    cb = colorbar; ylabel(cb, vars3D{k,3});
    title(vars3D{k,2});
    nz = ~isnan(C3D) & C3D>0;
xlim([min(x1(nz)) max(x1(nz))]);
ylim([min(y1(nz)) max(y1(nz))]);
zlim([min(z1(nz)) max(z1(nz))]);
clim([1e18 5e21])
end

try, colormap(turbo); catch, end



figure('Name','3D: aveSpyl (sputtering yield)','Position',[100 100 600 450]);

C3D = abs(aveSpyl);
C3D(~isfinite(C3D)) = NaN;

patch([x1 x2 x3]', [y1 y2 y3]', [z1 z2 z3]', C3D', ...
      'FaceColor','flat','EdgeColor','none','FaceAlpha',1,'EdgeAlpha',0.1);

axis equal; view(3); grid on;
xlabel('X'); ylabel('Y'); zlabel('Z');
cb = colorbar; ylabel(cb,'aveSpyl (dimensionless)');
title('3D: Average sputtering yield (aveSpyl)');
xlim([min(x1(nz)) max(x1(nz))]);
ylim([min(y1(nz)) max(y1(nz))]);
zlim([min(z1(nz)) max(z1(nz))]);
set(gca, 'ColorScale', 'log', 'FontSize', 16);

try, colormap(turbo); catch, end

%% --- Area-integrated totals (particles/s) ---
% Your per-area arrays are in particles/m^2/s, so multiplying by m^2 gives particles/s.

A = ElementArea;  % m^2 per triangle

grossEro_int  = nansum(grossPlasmaA  .* A);   % total gross erosion rate [particles/s]
grossDepo_int = nansum(grossDepoA .* A);   % total gross deposition rate [particles/s]
netEro_int    = nansum(netEroA    .* A);   % total net erosion rate [particles/s]

selfEro_int    = nansum(grossEroRawA    .* A);   % total net erosion rate [particles/s]


% (Optional) If you want net deposition as a positive number when deposition dominates:
netDepo_int   = -netEro_int;              % [particles/s], positive means net deposition

% Print results
fprintf('\n=== Area-integrated totals ===\n');
fprintf('Total area                 = %.6e m^2\n', nansum(A));
fprintf('Gross erosion (integrated) = %.6e particles/s\n', grossEro_int);
fprintf('Gross depo   (integrated)  = %.6e particles/s\n', grossDepo_int);
fprintf('Net erosion  (integrated)  = %.6e particles/s\n', netEro_int);
fprintf('Self erosion  (integrated)  = %.6e particles/s\n', selfEro_int);
fprintf('Net depo     (integrated)  = %.6e particles/s (positive => net deposition)\n', netDepo_int);