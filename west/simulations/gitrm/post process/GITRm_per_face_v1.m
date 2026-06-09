clc; clear; close all; format long

%% --- Load GITRm per-face data ---
% path1 = '../WEST_runs/gitrm-surface2.nc';
path1 = '../Guillaume_runs/1p1e17/WEST_runs/gitrm-surface2.nc';
% path1 = '../Guillaume_runs/7p3e18/WEST_runs/gitrm-surface2.nc';
ncid1 = netcdf.open(path1,'NOWRITE');

coordsID      = netcdf.inqVarID(ncid1,'Coordinates');
plasmaID      = netcdf.inqVarID(ncid1,'PlasmaErosion');
depoID        = netcdf.inqVarID(ncid1,'grossDeposition');
grossEroID    = netcdf.inqVarID(ncid1,'grossErosion');

coords        = netcdf.getVar(ncid1, coordsID);
grossPlasma   = netcdf.getVar(ncid1, plasmaID);
grossDepo     = netcdf.getVar(ncid1, depoID);
grossEroRaw   = squeeze(netcdf.getVar(ncid1, grossEroID));

netcdf.close(ncid1);

total_erosion = sum(grossPlasma);
totalDepostion = sum(grossDepo);
netErosion = sum(total_erosion-totalDepostion);
% Number of triangles
nTri = numel(grossPlasma);
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
% erosion_rate = 5.175509e+20 ; % [particles/s] - High Density WEST
erosion_rate = 2.021521e+19; % [particles/s] - Low Density WEST
nP           = 1e6;                    % primary particles
erosionPP    = erosion_rate / nP;

% Normalize per face by area (-> per-area rates/fluxes)
grossPlasmaA   = erosionPP .* grossPlasma   ./ ElementArea;
grossDepoA     = erosionPP .* grossDepo     ./ ElementArea;
grossEroRawA   = erosionPP .* grossEroRaw   ./ ElementArea;

% grossPlasmaA   = grossPlasma ;
% grossDepoA     = grossDepo   ;
% grossEroRawA   =  grossEroRaw ;

% If "gross erosion" should include plasma-driven component:
grossEroA = grossPlasmaA + grossEroRawA;   % per-area gross erosion
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

    axis equal
    xlabel('X'); ylabel('Y'); zlabel('Z');
    grid on; view(3);
    cb = colorbar; ylabel(cb, vars3D{k,3});
    title(vars3D{k,2});
% xlim([8.1 8.6]);
% ylim([-1 1]);
 % zlim([-1 1]);
set(gca, 'ColorScale', 'log', 'FontSize', 16);
clim([1e17 1e23])
end

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


return

%% ==========================================
% Geometry offsets (your values)
tor_x0 = 0; %1.49345169307177;
tor_y0 = 2.21; %0.011285617584263547;

% Cylindrical per vertex (match your convention: phi = atan2(X, Y))
R1 = sqrt(x1.^2 + y1.^2);  R2 = sqrt(x2.^2 + y2.^2);  R3 = sqrt(x3.^2 + y3.^2);
phi1 = atan2(x1, y1);       phi2 = atan2(x2, y2);       phi3 = atan2(x3, y3);
theta1 = atan2(z1 - tor_y0, R1 - tor_x0);
theta2 = atan2(z2 - tor_y0, R2 - tor_x0);
theta3 = atan2(z3 - tor_y0, R3 - tor_x0);

% Per-triangle unwrap of phi to avoid seam
phi_mat   = [phi1,   phi2,   phi3];    % sz x 3
theta_mat = [theta1, theta2, theta3];  % sz x 3
phi_unw   = phi_mat;
for k = 2:3
    d = phi_unw(:,k) - phi_unw(:,1);
    phi_unw(:,k) = phi_unw(:,k) - 2*pi*(d >  pi) + 2*pi*(d < -pi);
end

Phi   = phi_unw.';         % 3 x sz
Theta = theta_mat.';       % 3 x sz
Zzero = zeros(3, sz);

Cphth = log10(abs(grossEroA)).';   % 1 x sz
Cphth(~isfinite(Cphth)) = NaN;
%  ==========================================
% (Uses Phi, Theta, Zzero you built above)
varsCyl = { ...
    grossEroA,  '(a)  Gross Erosion ', '[Particles/m^2/s]'; ...
    grossDepoA, '(b)  Gross Deposition',    '[Particles/m^2/s]'; ...
    netEroA,    '(c) Net Erosion',   '[Particles/m^2/s]';
    netDepoA,    '(d) Net Deposition',   '[Particles/m^2/s]' ...
    };

figure('Name','(\phi,\theta) map (linear-scale, per area)','Position',[120 120 1200 450]);

for k = 1:4
    subplot(2,2,k)
    Cphth = ((varsCyl{k,1})).';
    Cphth(~isfinite(Cphth)) = NaN;

    patch(Phi, Theta, Zzero, Cphth, ...
          'FaceColor','flat', 'EdgeAlpha',0.15, 'FaceAlpha',1);

    view(2);  %axis equal tight;
    box on;
    % xlim([-0.484 0.484]);
    % ylim([-0.7032 -0.635])
    xlabel('\phi'); ylabel('\theta'); grid on
    cb2 = colorbar; ylabel(cb2, varsCyl{k,3});
    title(varsCyl{k,2});
end

try, colormap("jet"); catch, end


