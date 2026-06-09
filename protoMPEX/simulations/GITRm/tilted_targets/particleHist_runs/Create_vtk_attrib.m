clc
clear
close all
format long

%Gitrm History file
path1="45Degree/";
file1="gitrm-history.nc";
fid=path1+file1;
x = ncread(fid, 'x');
y = ncread(fid, 'y');
z = ncread(fid, 'z');
sp=ncread(fid, 'species');

nP=size(x,2);
nT=size(x,1);
x1=reshape(x,[1, nP*nT]);
y1=reshape(y,[1, nP*nT]);
z1=reshape(z,[1, nP*nT]);
sp1=reshape(sp,[1, nP*nT]);

%% Gitrm Initial file to know which particle has what mass
%mass_id=sp1(1,1:nP*(nT-1));
mass_id = sp(1:end-1, :);
mass_id = reshape(mass_id, [1, nP*(nT-1)]);

%% Create Vtk file for visualization
n_cells=nP*(nT-1);
cells=zeros(n_cells,3);
cell_type=zeros(n_cells,1);

cumul=0;
for i=1:nP
  for j=1:nT
    if j<=nT-1
      index=(i-1)*(nT-1)+j;
      cells(index,:)=[2 cumul cumul+1];
      cell_type(index)=3 ;
    end
    cumul=cumul+1 ;
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
v = VideoWriter(fullfile(path1, 'Fixed_Z_Horizontal.avi'), 'Motion JPEG AVI');
v.FrameRate = 5; 
open(v);

targetWidth = 560;
targetHeight = 420;

fig = figure('Color', 'w');
set(fig, 'Units', 'pixels', 'Position', [100, 100, targetWidth, targetHeight]);

ax = axes('Parent', fig);
grid(ax, 'off'); hold(ax, 'on');

% --- FIXING THE FRAME OF REFERENCE ---
% Calculate global limits once so the box never moves
x_lims = [min(x(:)) max(x(:))];
y_lims = [min(y(:)) max(y(:))];
z_lims = [min(z(:)) max(z(:))];

xlim(ax, x_lims);
ylim(ax, y_lims);
zlim(ax, z_lims);

% This ensures 1 unit in X is the same size as 1 unit in Z (no stretching)
daspect(ax, [1 1 1]); 

% LOCK the axes
axis(ax, 'manual');

% --- ROTATE VIEW: MAKE Z HORIZONTAL ---
% view(0, 0) looks at the X-Z plane. 
% In this view, Z will be the horizontal axis and X will be vertical.
view(ax, 100, 10); 

xlabel(ax, 'X'); ylabel(ax, 'Y'); zlabel(ax, 'Z');
%title(ax, 'Fixed Reference Frame: Z-axis Horizontal');

unique_sp = unique(sp(:));
colors = lines(max(length(unique_sp), 1));

frameStep = 1; 

for t = 1:frameStep:nT
    % Use a tag to identify paths so we can delete them without clearing the whole axes
    % This is faster and keeps the labels/grid perfectly still
    if t > 1
        delete(findall(ax, 'Tag', 'dynamic_path'));
    end
    
    for p = 1:nP
        sp_val = sp(t, p);
        c_idx = find(unique_sp == sp_val);
        if isempty(c_idx), c_idx = 1; end
        
        % Plotting history
        plot3(ax, x(1:t, p), y(1:t, p), z(1:t, p), ...
              'Color', colors(c_idx, :), 'LineWidth', 1, 'Tag', 'dynamic_path');
        
        % Plotting current head
        plot3(ax, x(t, p), y(t, p), z(t, p), 'o', ...
              'MarkerEdgeColor', 'k', 'MarkerFaceColor', colors(c_idx, :), ...
              'MarkerSize', 4, 'Tag', 'dynamic_path');
    end
    
    drawnow;
    
    % Frame Size Correction
    currentFrame = getframe(fig);
    fixedRGB = imresize(currentFrame.cdata, [targetHeight, targetWidth]);
    writeVideo(v, fixedRGB);
end

close(v);
fprintf('Video saved with fixed Z-horizontal perspective.\n');