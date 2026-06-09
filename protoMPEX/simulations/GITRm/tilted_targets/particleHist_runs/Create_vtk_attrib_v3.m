clc
clear
close all
format long

%% =========================
%  INPUT FILES
% =========================
path1    = "0Degree/";
histFile = path1 + "gitrm-history.nc";

% ---- Geometry STL (exported from ParaView) ----
% If you want to use your uploaded STL directly, set:
% geomStl = "/mnt/data/45_geom.stl";
geomStl  = path1 + "0_geom.stl";   % <=== change name/path if needed

%% =========================
%  READ HISTORY
% =========================
x  = ncread(histFile, 'x');
y  = ncread(histFile, 'y');
z  = ncread(histFile, 'z');
sp = ncread(histFile, 'species');

nP = size(x,2);
nT = size(x,1);

x1 = reshape(x,[1, nP*nT]);
y1 = reshape(y,[1, nP*nT]);
z1 = reshape(z,[1, nP*nT]);

%% =========================
%  MASS ID
% =========================
mass_id = sp(1:end-1, :);
mass_id = reshape(mass_id, [1, nP*(nT-1)]);

%% =========================
%  CREATE VTK (particle paths) [unchanged]
% =========================
n_cells   = nP*(nT-1);
cells     = zeros(n_cells,3);
cell_type = zeros(n_cells,1);

cumul = 0;
for i=1:nP
  for j=1:nT
    if j<=nT-1
      index=(i-1)*(nT-1)+j;
      cells(index,:)=[2 cumul cumul+1];
      cell_type(index)=3;
    end
    cumul=cumul+1;
  end
end

fileID = fopen(path1 + "XX1.vtk",'w');
fprintf(fileID,'%6s\n','# vtk DataFile Version 2.0');
fprintf(fileID,'%6s\n','particlepaths');
fprintf(fileID,'%4s\n','ASCII');
fprintf(fileID,'%4s\n','DATASET UNSTRUCTURED_GRID');
fprintf(fileID,'%4s %d %4s \n','POINTS', nP*nT, 'double');
fprintf(fileID,'%0.15f %0.15f %0.15f\n',[x1;y1;z1]);
fprintf(fileID,'%4s %d %d \n','CELLS', n_cells, 3*n_cells);
fprintf(fileID,'%d %d %d \n', cells');
fprintf(fileID,'%4s %d\n','CELL_TYPES', n_cells);
fprintf(fileID,'%d\n',cell_type);
fprintf(fileID,'%4s %d \n','CELL_DATA', nP*(nT-1));
fprintf(fileID,'SCALARS mass int 1\n');
fprintf(fileID,'LOOKUP_TABLE default\n');
fprintf(fileID,'%d\n',mass_id);
fclose(fileID);

%% ===========================
%   MOVIE (Mac-compatible MP4)
% ===========================
clc
v = VideoWriter(fullfile(path1, 'Fixed_Z_Horizontal_NIGHT.mp4'), 'MPEG-4'); % Mac compatible
v.FrameRate = 5;
v.Quality   = 95;   % 0–100
open(v);

targetWidth  = 560;
targetHeight = 420;

fig = figure('Color',[0.06 0.06 0.08]);
set(fig,'Units','pixels','Position',[100 100 targetWidth targetHeight]);

ax = axes('Parent',fig);
hold(ax,'on');

% ---- FIXED FRAME ----
x_lims = [min(x(:)) max(x(:))];
y_lims = [min(y(:)) max(y(:))];
z_lims = [min(z(:)) max(z(:))];

xlim(ax,x_lims);
ylim(ax,y_lims);
zlim(ax,z_lims);

daspect(ax,[1 1 1]);
axis(ax,'manual');

view(ax,100,10);

xlabel(ax,'X (m)','Color',[0.92 0.92 0.95]);
ylabel(ax,'Y (m)','Color',[0.92 0.92 0.95]);
zlabel(ax,'Z (m)','Color',[0.92 0.92 0.95]);

%% ===== NIGHT MODE AXES BOX =====
ax.Color = [0.06 0.06 0.08];
box(ax,'on');
grid(ax,'on');

ax.LineWidth = 3.0;

ax.XColor = [0.92 0.92 0.95];
ax.YColor = [0.92 0.92 0.95];
ax.ZColor = [0.92 0.92 0.95];

ax.GridColor = [1 1 1];
ax.GridAlpha = 0.12;

ax.TickDir = 'out';
ax.TickLength = [0.015 0.015];

%% ===== RENDERING =====
set(fig,'Renderer','opengl');
set(ax,'Projection','perspective');
axis(ax,'vis3d');

%% =========================
%  GEOMETRY PATCH FROM STL
%  (Opaque cutaway + slightly darker side surfaces)
% =========================
if isfile(geomStl)
    TR = stlread(geomStl);
    Vg = TR.Points;
    Fg = TR.ConnectivityList;

    % ---- Use geometry bounds for axes ----
    geomMin = min(Vg, [], 1);
    geomMax = max(Vg, [], 1);

    padFrac   = 0.05;
    geomRange = geomMax - geomMin;
    geomMin   = geomMin - padFrac * geomRange;
    geomMax   = geomMax + padFrac * geomRange;

    xlim(ax, [geomMin(1) geomMax(1)]);
    ylim(ax, [geomMin(2) geomMax(2)]);
    zlim(ax, [geomMin(3) geomMax(3)]);
    daspect(ax,[1 1 1]);
    axis(ax,'manual');
    axis(ax,'vis3d');

    % ---- Compute face normals ----
    v1 = Vg(Fg(:,2),:) - Vg(Fg(:,1),:);
    v2 = Vg(Fg(:,3),:) - Vg(Fg(:,1),:);
    fn = cross(v1, v2, 2);
    fn = fn ./ vecnorm(fn,2,2);

    % ---- Open the geometry towards camera (cutaway) ----
    openX = fn(:,1) > 0.6;
    openY = fn(:,2) > 0.6;
    keep  = ~(openX | openY);
    Fcut  = Fg(keep,:);

    % ---- Make side-ish faces slightly darker (based on normal) ----
    fnCut = fn(keep,:);
    nz    = abs(fnCut(:,3)); % vertical component; walls have nz ~ 0, floors have nz ~ 1
    isSide = nz < 0.35;      % "side walls" threshold (tune 0.25–0.45)

    F_side = Fcut(isSide,:);
    F_main = Fcut(~isSide,:);

    mainGray = [0.78 0.78 0.78];  % main body a bit brighter
    sideGray = [0.45 0.45 0.45];  % side walls darker

    % ---- MAIN (gray), two-sided ----
    hMain1 = patch(ax, ...
        'Faces', F_main, 'Vertices', Vg, ...
        'FaceColor', mainGray, ...
        'EdgeColor', 'none', ...
        'FaceAlpha', 1.0, ...
        'AmbientStrength', 0.20, ...
        'DiffuseStrength', 0.55, ...
        'SpecularStrength', 0.55, ...
        'SpecularExponent', 35, ...
        'BackFaceLighting', 'reverselit', ...
        'Tag', 'static_geom');

    hMain2 = patch(ax, ...
        'Faces', F_main(:,[1 3 2]), 'Vertices', Vg, ...
        'FaceColor', mainGray, ...
        'EdgeColor', 'none', ...
        'FaceAlpha', 1.0, ...
        'AmbientStrength', 0.20, ...
        'DiffuseStrength', 0.55, ...
        'SpecularStrength', 0.55, ...
        'SpecularExponent', 35, ...
        'BackFaceLighting', 'reverselit', ...
        'Tag', 'static_geom');

    % ---- SIDE WALLS (darker), two-sided ----
    if ~isempty(F_side)
        hSide1 = patch(ax, ...
            'Faces', F_side, 'Vertices', Vg, ...
            'FaceColor', sideGray, ...
            'EdgeColor', 'none', ...
            'FaceAlpha', 1.0, ...
            'AmbientStrength', 0.25, ...
            'DiffuseStrength', 0.60, ...
            'SpecularStrength', 0.35, ...
            'SpecularExponent', 25, ...
            'BackFaceLighting', 'reverselit', ...
            'Tag', 'static_geom');

        hSide2 = patch(ax, ...
            'Faces', F_side(:,[1 3 2]), 'Vertices', Vg, ...
            'FaceColor', sideGray, ...
            'EdgeColor', 'none', ...
            'FaceAlpha', 1.0, ...
            'AmbientStrength', 0.25, ...
            'DiffuseStrength', 0.60, ...
            'SpecularStrength', 0.35, ...
            'SpecularExponent', 25, ...
            'BackFaceLighting', 'reverselit', ...
            'Tag', 'static_geom');
    end

    % ---- Lighting ----
    lighting(ax,'gouraud');
    camlight(ax,'headlight');
    camlight(ax,'right');
    material(hMain1,'shiny');
    material(hMain2,'shiny');
else
    warning("STL geometry not found: %s", geomStl);
end

%% ===== PARTICLES (yellow with red edge) =====
tailColor = [1 1 0];    % yellow
headFace  = [1 1 0];    % yellow fill
headEdge  = [1 0 0];    % red edge
headSize  = 3.5;

%% ===== Reduce number of particles plotted =====
particleStep  = 10;               % plot every 10th particle
plotParticles = 1:particleStep:nP;

frameStep = 1;

for t = 1:frameStep:nT
    if t > 1
        delete(findall(ax,'Tag','dynamic_path'));
    end

    for p = plotParticles
        % Tail
        plot3(ax, x(1:t,p), y(1:t,p), z(1:t,p), ...
              'Color', tailColor, 'LineWidth', 1.2, 'Tag','dynamic_path');

        % Head
        plot3(ax, x(t,p), y(t,p), z(t,p), 'o', ...
              'MarkerFaceColor', headFace, ...
              'MarkerEdgeColor', headEdge, ...
              'LineWidth', 1.1, ...
              'MarkerSize', headSize, ...
              'Tag','dynamic_path');
    end

    drawnow;

    frame = getframe(fig); % no tight axes
    fixedRGB = imresize(frame.cdata,[targetHeight targetWidth]);
    writeVideo(v,fixedRGB);
end

close(v);
fprintf('Night-mode MP4 saved (Mac compatible) with STL geometry (cutaway) + darker side walls + reduced particles.\n');