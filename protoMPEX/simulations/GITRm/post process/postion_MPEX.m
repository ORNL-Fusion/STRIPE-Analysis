close all;
clear all;

data = readtable('../MPEX_runs/finalPos0.txt');

X = [data.x_m_, data.y_m_, data.z_m_];

% ---- density estimate using kNN ----
k = 20;   % neighbors (15–30 is a good range)
[idx, dist] = knnsearch(X, X, 'K', k+1);
rk = dist(:,end);                 % distance to k-th neighbor
density = k ./ ((4/3)*pi*rk.^3);  % local number density



% ---- plotting ----
figure('Color','w');
ax = axes; hold(ax,'on');

% dark background (optional but looks great)
ax.Color = [0.06 0.06 0.08];
ax.XColor = [0.9 0.9 0.9];
ax.YColor = [0.9 0.9 0.9];
ax.ZColor = [0.9 0.9 0.9];
grid on; box on;

% glow layer
scatter3(data.x_m_, data.y_m_, data.z_m_, ...
    18, density, 'filled', ...
    'MarkerFaceAlpha',0.08, 'MarkerEdgeAlpha',0);

% sharp layer
scatter3(data.x_m_, data.y_m_, data.z_m_, ...
    6, density, 'filled', ...
    'MarkerFaceAlpha',0.7, 'MarkerEdgeAlpha',0);

axis equal
xlabel('x [m]')
ylabel('y [m]')
zlabel('z [m]')

colormap(parula)
cb = colorbar;
cb.Label.String = 'Local particle density (arb.)';

set(ax,'Projection','perspective')
view(3)
lighting gouraud
camlight headlight
camlight right
material shiny
density = log10(density);
cb.Label.String = 'log_{10}(density)';