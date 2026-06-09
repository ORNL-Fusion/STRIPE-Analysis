close all;
clear all;

% Read table with default naming (so x_m_, y_m_, z_m_ exist)
data = readtable('../tilted_targets/test/0_degrees/MPEX_runs/finalPos0.txt');

x = data.x_m_;
y = data.y_m_;
z = data.z_m_;

% -------------------------
% 3D density via binning (no histcounts3 needed)
% -------------------------
nb = 60;

edgesX = linspace(min(x), max(x), nb+1);
edgesY = linspace(min(y), max(y), nb+1);
edgesZ = linspace(min(z), max(z), nb+1);

binX = discretize(x, edgesX);
binY = discretize(y, edgesY);
binZ = discretize(z, edgesZ);

inside = ~isnan(binX) & ~isnan(binY) & ~isnan(binZ);

N = accumarray([binX(inside), binY(inside), binZ(inside)], 1, [nb nb nb]);

density = nan(size(x));
linInd = sub2ind([nb nb nb], binX(inside), binY(inside), binZ(inside));
density(inside) = N(linInd);

% Optional (recommended): log scale for better contrast
useLog = true;
if useLog
    densityPlot = log10(density);
    cbLabel = 'log_{10}(binned particle count)';
else
    densityPlot = density;
    cbLabel = 'Binned particle count';
end

% -------------------------
% Fancy plot
% -------------------------
figure('Color','w');
ax = axes; hold(ax,'on');

ax.Color = [0.06 0.06 0.08];
ax.GridColor = [1 1 1];
ax.GridAlpha = 0.12;
ax.XColor = [0.9 0.9 0.92];
ax.YColor = [0.9 0.9 0.92];
ax.ZColor = [0.9 0.9 0.92];
box on; grid on;

% Glow layer
scatter3(x(inside), y(inside), z(inside), ...
    18, densityPlot(inside), 'filled', ...
    'MarkerFaceAlpha', 0.08, 'MarkerEdgeAlpha', 0);

% Crisp layer
scatter3(x(inside), y(inside), z(inside), ...
    6, densityPlot(inside), 'filled', ...
    'MarkerFaceAlpha', 0.7, 'MarkerEdgeAlpha', 0);

axis equal
xlabel('x [m]','FontSize',12,'Color',[0.92 0.92 0.95])
ylabel('y [m]','FontSize',12,'Color',[0.92 0.92 0.95])
zlabel('z [m]','FontSize',12,'Color',[0.92 0.92 0.95])

colormap(parula)
cb = colorbar;
cb.Color = [0.92 0.92 0.95];
cb.Label.String = cbLabel;

valid = isfinite(densityPlot);
if any(valid)
    clim = prctile(densityPlot(valid), [2 98]);  % robust limits
    if clim(1) == clim(2)
        clim = [min(densityPlot(valid)) max(densityPlot(valid))];
    end
    caxis(clim)
end

view(3)
set(ax,'Projection','perspective')
lighting gouraud
camlight headlight
camlight right
material shiny

title('Particle Final Positions (colored by density)', ...
      'FontWeight','normal','Color',[0.95 0.95 0.98]);