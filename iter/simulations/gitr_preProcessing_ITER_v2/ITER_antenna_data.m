%% --- Constants ---
ME = 9.10938356e-31;
MI = 1.6737236e-27;
Q = 1.60217662e-19;
EPS0 = 8.854187e-12;

amu = 18;  % Atomic mass unit of the species

%% --- Read profiles from SOLEDGE ---
% Grid
R = ncread('profilesITER.nc', 'x');
z = ncread('profilesITER.nc', 'z');

% Ensure R and z are meshgrids for proper interpolation
[R_grid, Z_grid] = meshgrid(R, z);

% B-field
bz = ncread('profilesITER.nc', 'bz');
bt = ncread('profilesITER.nc', 'bt');
br = ncread('profilesITER.nc', 'br');

% Densities
ne = ncread('profilesITER.nc', 'ne');
ni = ncread('profilesITER.nc', 'ni');

% Temperatures
te = ncread('profilesITER.nc', 'te');
ti = ncread('profilesITER.nc', 'ti');

% Velocity
vt = ncread('profilesITER.nc', 'vt');
vr = ncread('profilesITER.nc', 'vr');
vz = ncread('profilesITER.nc', 'vz');

% Ensure vz has the same shape as vr and vt
if numel(vz) ~= numel(vr)
    vz = reshape(vz, size(vr)); % Reshape if necessary
end

% Compute centroids
r_centroid = sqrt(centroid(:,1).^2 + centroid(:,2).^2);

% Plot ni profile
figure; imagesc(R, z, ni);
set(gca, 'YDir', 'normal', 'FontName', 'times', 'FontSize', 18);
ylabel('$r$ [m]', 'Interpreter', 'Latex', 'FontSize', 18);
xlabel('$z$ [m]', 'Interpreter', 'Latex', 'FontSize', 18);
title('Input ni');
colorbar;
hold on;
scatter(r_centroid, centroid(:,3), 'r');

% Plot vz profile
figure; imagesc(R, z, vz);
set(gca, 'YDir', 'normal', 'FontName', 'times', 'FontSize', 18);
ylabel('$r$ [m]', 'Interpreter', 'Latex', 'FontSize', 18);
xlabel('$z$ [m]', 'Interpreter', 'Latex', 'FontSize', 18);
title('Input vz');
colorbar;
hold on;
scatter(r_centroid, centroid(:,3), 'r');

%% --- Interpolations onto the geometry ---
interp_method = 'natural';  % Use spline interpolation for smoother results

ne_surf = interp2(R_grid, Z_grid, ne, r_centroid, centroid(:,3), interp_method);
ni_surf = interp2(R_grid, Z_grid, ni, r_centroid, centroid(:,3), interp_method);
te_surf = interp2(R_grid, Z_grid, te, r_centroid, centroid(:,3), interp_method);
ti_surf = interp2(R_grid, Z_grid, ti, r_centroid, centroid(:,3), interp_method);

br_surf = interp2(R_grid, Z_grid, br, r_centroid, centroid(:,3), interp_method);
bt_surf = interp2(R_grid, Z_grid, bt, r_centroid, centroid(:,3), interp_method);
bz_surf = interp2(R_grid, Z_grid, bz', r_centroid, centroid(:,3), interp_method);

vr_surf = interp2(R_grid, Z_grid, vr, r_centroid, centroid(:,3), interp_method);
vt_surf = interp2(R_grid, Z_grid, vt, r_centroid, centroid(:,3), interp_method);
vz_surf = interp2(R_grid, Z_grid, vz, r_centroid, centroid(:,3), interp_method);

% Replace NaN values with zero
fields = {'ne_surf', 'ni_surf', 'te_surf', 'ti_surf', ...
          'vr_surf', 'vz_surf', 'vt_surf', 'br_surf', 'bz_surf', 'bt_surf'};
for i = 1:length(fields)
    eval([fields{i} '(isnan(' fields{i} ')) = 0;']);
end

%% --- Magnetic field transformation ---
phi_centroid = atan2(centroid(:,2), centroid(:,1));

bx = br_surf .* cos(phi_centroid) - bt_surf .* sin(phi_centroid);
by = br_surf .* sin(phi_centroid) + bt_surf .* cos(phi_centroid);
bz = bz_surf;

b_mag = sqrt(bx.^2 + by.^2 + bz.^2);
ubx = bx ./ b_mag;
uby = by ./ b_mag;
ubz = bz ./ b_mag;

%% --- Velocity transformation ---
vx = vr_surf .* cos(phi_centroid) - vt_surf .* sin(phi_centroid);
vy = vr_surf .* sin(phi_centroid) + vt_surf .* cos(phi_centroid);
vz = vz_surf;

v_mag = sqrt(vx.^2 + vy.^2 + vz.^2);
uvx = vx ./ v_mag;
uvy = vy ./ v_mag;
uvz = vz ./ v_mag; % Fixed incorrect usage of vy

% Oxygen flux
o8plus_flux_surf = ni_surf .* v_mag;

%% --- Normal vectors and theta calculation ---
norm_vec_mag = sqrt(norm_vec(:,1).^2 + norm_vec(:,2).^2 + norm_vec(:,3).^2);
unorm_vec = norm_vec ./ norm_vec_mag;

theta = acos(unorm_vec(:,1) .* ubx + unorm_vec(:,2) .* uby + unorm_vec(:,3) .* ubz);
theta(isnan(theta)) = 0;

% Ensure theta is within [0, pi/2]
ii = find(theta > pi/2);
% theta(ii) = abs(theta(ii) - pi);
theta(ii) = (theta(ii) - pi);

%% --- Save Outputs ---
writematrix(o8plus_flux_surf, 'o8plus_flux_surf.csv');

%% --- Visualization ---
% Normal and B-field vectors
figure; 
quiver3(centroid(:,1), centroid(:,2), centroid(:,3), norm_vec(:,1), norm_vec(:,2), norm_vec(:,3)); 
hold on;
quiver3(centroid(:,1), centroid(:,2), centroid(:,3), bx, by, bz);

% Histogram of theta
figure; histogram(theta);
% 
% figure; 
% patch(transpose(X), transpose(Y), transpose(Z), v_mag, ...
%     'FaceAlpha', 1, 'EdgeAlpha', 0.3);
% title('Flow on surface', 'FontSize', 18);
% colorbar('eastoutside');
% xlabel('X', 'FontSize', 18);
% ylabel('Y', 'FontSize', 18);
% zlabel('Z', 'FontSize', 18);
% set(gca, 'FontSize', 18);
