clc
clear
close all
format long

%Gitrm History file
path1="45Degree/";
file1="gitrm-history.nc";
fid=path1+file1;

x  = ncread(fid, 'x');
y  = ncread(fid, 'y');
z  = ncread(fid, 'z');
sp = ncread(fid, 'species');

nP = size(x,2);
nT = size(x,1);

x1 = reshape(x,[1, nP*nT]);
y1 = reshape(y,[1, nP*nT]);
z1 = reshape(z,[1, nP*nT]);

%% Gitrm Initial file to know which particle has what mass
mass_id = sp(1:end-1, :);
mass_id = reshape(mass_id, [1, nP*(nT-1)]);

%% Create Vtk file for visualization
n_cells = nP*(nT-1);
cells = zeros(n_cells,3);
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

%%
clc
v = VideoWriter(fullfile(path1, 'Fixed_Z_Horizontal_NIGHT.avi'), 'Motion JPEG AVI');
v.FrameRate = 5; 
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

%% ===== PARTICLES (NIGHT MODE COLORS) =====
tailColor = [1 1 0];    % bright yellow
headFace  = [1 1 0];    % yellow fill
headEdge  = [1 0 0];    % red edge
headSize  = 7;

frameStep = 1;

for t = 1:frameStep:nT
    if t > 1
        delete(findall(ax,'Tag','dynamic_path'));
    end
    
    for p = 1:nP
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
    
    frame = getframe(fig);
    fixedRGB = imresize(frame.cdata,[targetHeight targetWidth]);
    writeVideo(v,fixedRGB);
end

close(v);
fprintf('Night-mode video saved.\n');