clc
clear 
close all
format long


%% GITRm per face
path1='../WEST_runs_test/gitrm-surface2.nc';
ncid1 = netcdf.open(path1,'NOWRITE');
coords=netcdf.getVar(ncid1,2);
grossPlasma=netcdf.getVar(ncid1,3);
grossDepo=netcdf.getVar(ncid1,4);
grossEro=netcdf.getVar(ncid1,5);
% GE_r1=netcdf.getVar(ncid1,4);
sz=size(grossEro,1);

x1=zeros(sz,1);
x2=zeros(sz,1);
x3=zeros(sz,1);
y1=zeros(sz,1);
y2=zeros(sz,1);
y3=zeros(sz,1);
z1=zeros(sz,1);
z2=zeros(sz,1);
z3=zeros(sz,1);

for i=1:sz
x1(i)=coords((i-1)*9+1);
y1(i)=coords((i-1)*9+2);
z1(i)=coords((i-1)*9+3);
x2(i)=coords((i-1)*9+4);
y2(i)=coords((i-1)*9+5);
z2(i)=coords((i-1)*9+6);
x3(i)=coords((i-1)*9+7);
y3(i)=coords((i-1)*9+8);
z3(i)=coords((i-1)*9+9);

end

cent_x=(x1+x2+x3)/3;
cent_y=(y1+y2+y3)/3;
cent_z=(z1+z2+z3)/3;

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
ElementArea = 0.5 * sqrt(sum(cp.^2, 2));  % area is sz x 1
% Erosion pre-factors
%% Erosion pre-factors
erosion_rate=5.9E19;
nP=1e6;
erosionPP=erosion_rate/nP;

grossPlasma = erosionPP.*grossPlasma./ElementArea;
grossParticle = erosionPP.*(grossEro)./ElementArea;
grossDepo = erosionPP.*grossDepo./ElementArea;

% grossEro=grossPlasma+grossEro;
grossEro=grossPlasma+grossParticle;

netEro=grossEro-grossDepo;

figure;
patch( ...
    [x1(:) x2(:) x3(:)]', ...   % X coordinates (3×Nfaces)
    [y1(:) y2(:) y3(:)]', ...   % Y coordinates (3×Nfaces)
    [z1(:) z2(:) z3(:)]', ...   % Z coordinates (3×Nfaces)
    grossEro, ...      % 1×Nfaces color data
    'FaceColor', 'flat', ...    % flat shading (per face)
    'EdgeColor', 'none');       % remove edges for cleaner look

axis equal
xlabel('X'); ylabel('Y'); zlabel('Z');
grid on; view(3);
colorbar
xlim([8.1 8.6]);
ylim([-1 1]);
zlim([-1 2]);
set(gca, 'ColorScale', 'log', 'FontSize', 16);

figure;
patch( ...
    [x1(:) x2(:) x3(:)]', ...   % X coordinates (3×Nfaces)
    [y1(:) y2(:) y3(:)]', ...   % Y coordinates (3×Nfaces)
    [z1(:) z2(:) z3(:)]', ...   % Z coordinates (3×Nfaces)
    abs(grossDepo), ...      % 1×Nfaces color data
    'FaceColor', 'flat', ...    % flat shading (per face)
    'EdgeColor', 'none');       % remove edges for cleaner look

axis equal
xlabel('X'); ylabel('Y'); zlabel('Z');
grid on; view(3);
colorbar
% xlim([8.1 8.6]);
% ylim([-1 1]);
% zlim([-1 2]);
set(gca, 'ColorScale', 'log', 'FontSize', 16);

figure;
patch( ...
    [x1(:) x2(:) x3(:)]', ...   % X coordinates (3×Nfaces)
    [y1(:) y2(:) y3(:)]', ...   % Y coordinates (3×Nfaces)
    [z1(:) z2(:) z3(:)]', ...   % Z coordinates (3×Nfaces)
     netEro, ...      % 1×Nfaces color data
    'FaceColor', 'flat', ...    % flat shading (per face)
    'EdgeColor', 'none');       % remove edges for cleaner look

axis equal
xlabel('X'); ylabel('Y'); zlabel('Z');
grid on; view(3);
colorbar
% xlim([8.1 8.6]);
% ylim([-1 1]);
% zlim([-1 2]);
set(gca, 'ColorScale', 'log', 'FontSize', 16);