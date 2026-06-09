% Load the COMSOL data
comsol_data = readmatrix('../comsol/feb_high1800_VDC.csv');
xx = comsol_data(:,1);
yy = comsol_data(:,2);
zz = comsol_data(:,3);
ne_surf = comsol_data(:,4);
emag = comsol_data(:,5);

% Coordinate transformation
R1 = sqrt(xx.^2 + yy.^2);
tor_x0 = 1.49345169307177;
tor_y0 = 0.011285617584263547;
phi = atan2(xx, yy);
theta = atan2(zz - tor_y0, R1 - tor_x0);

% Remove duplicates explicitly
coords = [phi, theta];
[coords_unique, unique_idx] = unique(coords, 'rows', 'stable');
phi_u = coords_unique(:,1);
theta = coords_unique(:,2);
ne_surf_u = ne_surf(unique_idx);
emag_u = emag(unique_idx);

% Delaunay triangulation
tri = delaunay(coords_unique(:,1), coords_unique(:,2));

% Plot electron density using trisurf
figure;
trisurf(tri, coords_unique(:,1), coords_unique(:,2), zeros(size(coords_unique,1),1), ne_surf(unique_idx),...
    'EdgeColor','none','FaceColor','interp');
view(2); % Top-down view
colorbar;
colormap('parula');
xlabel('\phi (rad)');
ylabel('\theta (rad)');
title('Surface Electron Density (n_e)');
axis tight;

% Plot electric field magnitude using trisurf
figure;
trisurf(tri, coords_unique(:,1), coords_unique(:,2), zeros(size(coords_unique,1),1), emag(unique_idx),...
    'EdgeColor','none','FaceColor','interp');
view(2); % Top-down view
colorbar;
colormap('parula');
xlabel('\phi (rad)');
ylabel('\theta (rad)');
title('Electric Field Magnitude (emag)');
axis tight;