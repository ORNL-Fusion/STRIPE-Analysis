clc
clear
close all

%% Files
geomFile = 'gitrGeometryPointPlane3d.cfg';
histFile = "../code_tests_t10_history/gitrm-history.nc";

%% Load particle history
x = ncread(histFile,'x');
y = ncread(histFile,'y');
z = ncread(histFile,'z');

nT = size(x,1);
nP = size(x,2);

%% Load antenna geometry
fid = fopen(geomFile);
for i = 1:20
    tline = fgetl(fid);
    if i>2 && ischar(tline)
        evalc(tline);
    end
end
fclose(fid);

X = [x1(:) x2(:) x3(:)];
Y = [y1(:) y2(:) y3(:)];
Z = [z1(:) z2(:) z3(:)];

%% Create figure
fig = figure('Color','w');
ax = axes;
hold on

%% Antenna surface
patch('XData',X',...
      'YData',Y',...
      'ZData',Z',...
      'FaceColor',[0.7 0.7 0.7],...
      'FaceAlpha',0.25,...
      'EdgeColor','none');

camlight
lighting gouraud

%% Fixed domain
xlim([8.1 8.6]);
ylim([-1 1]);
zlim([-1 2]);

axis equal
view(3)

xlabel('X')
ylabel('Y')
zlabel('Z')

%% Particle scatter (red shiny points)
h = scatter3(nan,nan,nan,...
    15,...                  % particle size
    'filled',...
    'MarkerFaceColor','r',...
    'MarkerEdgeColor','k');

%% Video
v = VideoWriter('particle_movie.avi');
v.FrameRate = 10;
open(v);

%% Time loop
for t = 1:nT

    set(h,...
        'XData',x(t,:),...
        'YData',y(t,:),...
        'ZData',z(t,:));

    drawnow

    frame = getframe(fig);
    writeVideo(v,frame);

end

close(v)

disp('Movie saved.')