clc; clear; close all; format long

%% --- Load GITRm per-face data ---
path1 = '../ITER_antenna_runs/gitrm-surface2.nc';
ncid1 = netcdf.open(path1,'NOWRITE');
coords      = netcdf.getVar(ncid1, 2);   % flattened [9*sz x 1]: (x1 y1 z1 x2 y2 z2 x3 y3 z3) per tri
grossPlasma = netcdf.getVar(ncid1, 3);   % sz x 1
grossDepo   = netcdf.getVar(ncid1, 4);   % sz x 1
grossEroRaw = netcdf.getVar(ncid1, 5);   % sz x 1   (raw gross erosion counts)
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
erosion_rate= 2.2e+18; % # 200882
nP           = 5e6;                    % primary particles
erosionPP    = erosion_rate / nP;

% Normalize per face by area (-> per-area rates/fluxes)
grossPlasmaA   = erosionPP .* grossPlasma   ./ ElementArea;
grossDepoA     = erosionPP .* grossDepo     ./ ElementArea;
grossEroRawA   = erosionPP .* grossEroRaw   ./ ElementArea;

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
    nz = ~isnan(C3D) & C3D>0;
xlim([min(x1(nz)) max(x1(nz))]);
ylim([min(y1(nz)) max(y1(nz))]);
zlim([min(z1(nz)) max(z1(nz))]);
end

try, colormap(turbo); catch, end


%% ==========================================
% Geometry offsets (your values)
tor_x0 = 8.34; %1.49345169307177;
tor_y0 = 0.93; %0.011285617584263547;

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


