clc
clear
close all
format long

%% GITRm per face
path1 = '../tilted_targets/test/0_degrees/MPEX_runs_v1/gitrm-surface2.nc';

ncid1 = netcdf.open(path1,'NOWRITE');

% --- Read coordinates and erosion by variable name ---
coordsVarID = netcdf.inqVarID(ncid1,'Coordinates');
coords      = netcdf.getVar(ncid1,coordsVarID);

erosVarID   = netcdf.inqVarID(ncid1,'PlasmaErosion');
PlasmaErosion = netcdf.getVar(ncid1,erosVarID);  % [nMeshFaces x nSpecies]

% Ensure double for math
PlasmaErosion = double(PlasmaErosion);

erosVarID   = netcdf.inqVarID(ncid1,'grossErosion');
GrossErosion = netcdf.getVar(ncid1,erosVarID);  % [nMeshFaces x nSpecies]

% Ensure double for math
PlasmaErosion = double(PlasmaErosion);
GrossErosion = double(GrossErosion);

erosVarID   = netcdf.inqVarID(ncid1,'grossDeposition');
GrossDeposition = netcdf.getVar(ncid1,erosVarID);  % [nMeshFaces x nSpecies]
GrossDeposition= double(GrossDeposition);

% --- Use FIRST COLUMN as erosion ---
% PlasmaErosion(:,1) = erosion for species 1 (what you asked for)
grossEro_Ta = PlasmaErosion(:,1) + GrossErosion(:,1);
grossEro_W = PlasmaErosion(:,2) + GrossErosion(:,2);
grossDepo_Ta = GrossDeposition(:,1);
% grossEro_W = PlasmaErosion(:,2) + GrossErosion(:,2);
% Optional: use second column as deposition (if present),
% otherwise set to zero so code still runs.
% if size(PlasmaErosion,2) >= 2
%     grossDepo = PlasmaErosion(:,2);
% else
%     grossDepo = zeros(size(grossEro));
% end

sz = size(grossEro_Ta,1);

% ---------------- Geometry (unchanged) ----------------
x1 = zeros(sz,1); x2 = x1; x3 = x1;
y1 = x1; y2 = x1; y3 = x1;
z1 = x1; z2 = x1; z3 = x1;

for i = 1:sz
    base = (i-1)*9;
    x1(i) = coords(base+1);
    y1(i) = coords(base+2);
    z1(i) = coords(base+3);
    x2(i) = coords(base+4);
    y2(i) = coords(base+5);
    z2(i) = coords(base+6);
    x3(i) = coords(base+7);
    y3(i) = coords(base+8);
    z3(i) = coords(base+9);
end

cent_x = (x1 + x2 + x3)/3;
cent_y = (y1 + y2 + y3)/3;
cent_z = (z1 + z2 + z3)/3;

% Build vertex coordinate vectors
P1 = [x1(:), y1(:), z1(:)];
P2 = [x2(:), y2(:), z2(:)];
P3 = [x3(:), y3(:), z3(:)];

% Two edge vectors for each triangle
v1 = P2 - P1;
v2 = P3 - P1;

% Cross product for each triangle
cp = cross(v1, v2, 2);   % cross along rows

% Triangle area = half the magnitude of the cross product
ElementArea = 0.5 * sqrt(sum(cp.^2, 2));  % sz x 1

%% Erosion pre-factors
erosion_rate = 1.165799127631290e+16;
nP           = 5e6;
erosionPP    = erosion_rate / nP;

% Convert to [atoms / m^2 / s] (or whatever your units are)
grossEro_Ta  = erosionPP .* grossEro_W  ./ ElementArea;
grossDepo_Ta = erosionPP .* grossDepo_Ta ./ ElementArea;

% Net erosion = erosion - deposition
% netEro = grossEro - grossDepo;

%% Plots
%% Plots
figure;
patch( ...
    [x1(:) x2(:) x3(:)]', ...
    [y1(:) y2(:) y3(:)]', ...
    [z1(:) z2(:) z3(:)]', ...
    grossDepo_Ta, ...                % deposition values
    'FaceColor','flat', ...
    'EdgeColor','none');

axis equal
xlabel('X'); ylabel('Y'); zlabel('Z');
grid on; view(3);

cb = colorbar;
ylabel(cb, 'Particles m^{-2} s^{-1}', 'FontSize', 18);

title('85^o');
% zlim([0.1 0.2])

% figure;
% patch( ...
%     [x1(:) x2(:) x3(:)]', ...
%     [y1(:) y2(:) y3(:)]', ...
%     [z1(:) z2(:) z3(:)]', ...
%     netEro, ...                   % net erosion
%     'FaceColor','flat', ...
%     'EdgeColor','none');
% 
% axis equal
% xlabel('X'); ylabel('Y'); zlabel('Z');
% grid on; view(3);
% colorbar
% title('Net erosion (col 1 - col 2 of PlasmaErosion)')