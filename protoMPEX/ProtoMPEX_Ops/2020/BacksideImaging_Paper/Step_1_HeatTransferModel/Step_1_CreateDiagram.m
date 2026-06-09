% Step 1: Create figures for heat transfer model
 
clear all
close all
clc

saveFig = 1;

% Plate dimensions:
Lx = 10/100;
Ly = 10/100;
Lz = 3/100;

% Outer square dimensions:
x_edge = 1/100;
y_edge = 1/100;
z_edge = 1/100;
Lx_outer = Lx + x_edge;
Ly_outer = Ly + y_edge;
Lz_outer = Lz + z_edge;

% Create plate polygon:
plate.x = [-1,+1,+1,-1,-1]*Lx/2;
plate.y = [-1,-1,+1,+1,-1]*Ly/2;
plate.z = [-1,-1,+1,+1,-1]*Lz/2;
plate.color = [0.80,0.80,0.80];

% Heat flux profile
heatFlux.sideView.y = linspace(-Ly,+Ly)/2;
heatFlux.sideView.z = gaussian(heatFlux.sideView.y,0,Ly/7);
heatFlux.sideView.z = heatFlux.sideView.z/max(heatFlux.sideView.z);

% Create hatching:
% Hatching dimensions:
m = Ly/Lx;
x = linspace(-1,1,1e3)*Lx_outer;
y = linspace(-1,1,1e3)*Ly_outer;
z = linspace(-Lz,Lz_outer,1e3);

nx = 30;
ny = 30;
nz = 30;

% "x" intercepts:
dx = Lx_outer/nx;
x0 = 0:dx:(Lx_outer-dx);

% "y" intercepts:
dy = Ly_outer/ny;
y0 = 0:dy:(Ly_outer-dy);

% Line functions:
yy = @(xx,ii,jj) y0(ii) + m*(xx - x0(jj));


%% Font view figure:
% =========================================================================

close all

figure('color','w'); 
h{1}.figure = gcf;
h{1}.axes(1) = gca;

hold on; 

% Hatching:
kk = 1;
for ii = 1:numel(x0)
    for jj = 1:numel(y0)
        ff = yy(x,ii,jj);
        rng = find(abs(ff) <= Ly_outer/2 & abs(x) <= Lx_outer/2);        
        h{1}.hatch(kk) = plot(x(rng),ff(rng),'k');
        kk = kk + 1;
    end
end
set(h{1}.hatch,'LineWidth',1)

% Plate:
h{1}.plate = fill(plate.x,plate.y,plate.color);
set(h{1}.plate,'EdgeColor',plate.color)

% Coordinate system:
% Text arrows fields:
fields.HorizontalAlignment = 'center';
fields.Color = 'k';
fields.FontSize = 11;
fields.Interpreter = 'Latex';
fields.HeadLength = 7;
fields.HeadWidth = 7;
fields.HeadStyle = 'vback2';
fields.LineWidth = 2;

% "x" coordinate:
fields.String = '' ;
x =      [0,1/100];
y =      [0,0    ];
h{1}.coordinate(1) = myTextArrow(gca,x,y,fields);
h{1}.coordinateText(1) = text(0.02,1.3e-4,'$\hat{x}_*$');

% "y" coordinate:
fields.String = '' ;
y =      [0,1.25/100];
x =      [0,0    ];
h{1}.coordinate(2) = myTextArrow(gca,x,y,fields);
h{1}.coordinateText(2) = text(-0.0001,0.022,'$\hat{y}_*$');

% Formatting coordinate:
set(h{1}.coordinateText,'interpreter','latex','FontSize',15)

% Origin:
h{1}.origin = plot(0,0,'ko','MarkerFaceColor','k');

% Get scale factor:
h{1}.plotBoxAspectRatio = get(h{1}.axes,'PlotBoxAspectRatio');
scaleFactor.x = h{1}.plotBoxAspectRatio(2);
scaleFactor.y = h{1}.plotBoxAspectRatio(1);

% Formatting:
axis image
xlim([-1,1]*Lx_outer*0.7)
ylim([-1,1]*Ly_outer*0.7)
set(h{1}.axes(1),'YTick',[],'XTick',[],'XColor','w','YColor','w')

% Dimensions 1:
% -------------------------------------------------------------------------
kk = 1;
% thin line 1
nn = 1;
dim{kk}.thinline(nn).x1 = max(plate.x) + 0.2/100;
dim{kk}.thinline(nn).x2 = dim{kk}.thinline(nn).x1 + 1.5/100;
dim{kk}.thinline(nn).y1 = max(plate.y);
dim{kk}.thinline(nn).y2 = max(plate.y);
x1 = dim{kk}.thinline(nn).x1;
x2 = dim{kk}.thinline(nn).x2;
y1 = dim{kk}.thinline(nn).y1;
y2 = dim{kk}.thinline(nn).y2;
h{1}.dim{kk}.thinline(nn) = line([x1,x2],[y1,y2],'color','k');

% thin line 2
nn = 2;
dim{kk}.thinline(nn).x1 = max(plate.x) + 0.2/100;
dim{kk}.thinline(nn).x2 = dim{kk}.thinline(nn).x1 + 1.5/100;
dim{kk}.thinline(nn).y1 = min(plate.y);
dim{kk}.thinline(nn).y2 = min(plate.y);
x1 = dim{kk}.thinline(nn).x1;
x2 = dim{kk}.thinline(nn).x2;
y1 = dim{kk}.thinline(nn).y1;
y2 = dim{kk}.thinline(nn).y2;
h{1}.dim{kk}.thinline(nn) = line([x1,x2],[y1,y2],'color','k');

% Double arrow:
dim{kk}.doubleArrow.x1 = dim{kk}.thinline(nn).x2 - 0.2/100;
dim{kk}.doubleArrow.x2 = dim{kk}.doubleArrow.x1;
dim{kk}.doubleArrow.y1 = max(plate.y);
dim{kk}.doubleArrow.y2 = min(plate.y);

x1 = dim{kk}.doubleArrow.x1;
x2 = dim{kk}.doubleArrow.x2;
y1 = dim{kk}.doubleArrow.y1;
y2 = dim{kk}.doubleArrow.y2;
[figx,figy] = dsxy2figxy(gca,[x1 x2]*scaleFactor.x,[y1 y2]*scaleFactor.y);
h{1}.dim{kk}.darrow = annotation('doublearrow',figx,figy,'Color','k');
h{1}.dim{kk}.darrow.HeadStyle = 'vback2';

h{1}.dim{kk}.text = text(0.07,2e-3,'$L_y$');
set(h{1}.dim{kk}.text,'interpreter','latex','FontSize',15)

% Dimensions 2:
% -------------------------------------------------------------------------
fields.LineWidth = 2;
kk = 2;
% thin line 1
nn = 1;
dim{kk}.thinline(nn).y1 = min(plate.y) - 0.2/100;
dim{kk}.thinline(nn).y2 = dim{kk}.thinline(nn).y1 - 1.5/100;
dim{kk}.thinline(nn).x1 = min(plate.x);
dim{kk}.thinline(nn).x2 = min(plate.x);
x1 = dim{kk}.thinline(nn).x1;
x2 = dim{kk}.thinline(nn).x2;
y1 = dim{kk}.thinline(nn).y1;
y2 = dim{kk}.thinline(nn).y2;
h{1}.dim{kk}.thinline(nn) = line([x1,x2],[y1,y2],'color','k');

% thin line 2
nn = 2;
dim{kk}.thinline(nn).y1 = min(plate.y) - 0.2/100;
dim{kk}.thinline(nn).y2 = dim{kk}.thinline(nn).y1 - 1.5/100;
dim{kk}.thinline(nn).x1 = max(plate.x);
dim{kk}.thinline(nn).x2 = max(plate.x);
x1 = dim{kk}.thinline(nn).x1;
x2 = dim{kk}.thinline(nn).x2;
y1 = dim{kk}.thinline(nn).y1;
y2 = dim{kk}.thinline(nn).y2;
h{1}.dim{kk}.thinline(nn) = line([x1,x2],[y1,y2],'color','k');

% Double arrow:
dim{kk}.doubleArrow.y1 = dim{kk}.thinline(nn).y2 + 0.2/100;
dim{kk}.doubleArrow.y2 = dim{kk}.doubleArrow.y1;
dim{kk}.doubleArrow.x1 = max(plate.x);
dim{kk}.doubleArrow.x2 = min(plate.x);

x1 = dim{kk}.doubleArrow.x1;
x2 = dim{kk}.doubleArrow.x2;
y1 = dim{kk}.doubleArrow.y1;
y2 = dim{kk}.doubleArrow.y2;
[figx,figy] = dsxy2figxy(gca,[x1 x2]*scaleFactor.x,[y1 y2]*scaleFactor.y);
h{1}.dim{kk}.darrow = annotation('doublearrow',figx,figy,'Color','k');
h{1}.dim{kk}.darrow.HeadStyle = 'vback2';

h{1}.dim{kk}.text{kk} = text(2e-3,-0.073,'$L_x$');
set(h{1}.dim{kk}.text{kk},'interpreter','latex','FontSize',15)

if 0
    % Circles:
    p0.x = 0;
    p0.y = 0;
    rad = linspace(0.1,1,4)*max(plate.x);
    for cc = 1:numel(rad)
        circle{cc} = DrawCircle(p0,rad(cc),100);
        h{1}.circle{cc} = plot(circle{cc}.x,circle{cc}.y,'r','lineWidth',2);
    end
end

% handles to use on nex figure:
hgroup{1} = [h{1}.plate,h{1}.hatch];
hannot{1} = [h{1}.dim{1}.thinline,h{1}.dim{1}.darrow];

% Save figure:
% =========================================================================
if saveFig
    figureName = 'Step_1_FrontViewPlate';
    % PDF figure:
    exportgraphics(gcf,[figureName,'.pdf'],'Resolution',600); %,'ContentType', 'vector') 

    % TIFF figure:
    exportgraphics(gcf,[figureName,'.tiff'],'Resolution',600) 
end
%% Side view:
% =========================================================================
% Copy the previous figure
figure(h{1}.figure);
h{2}.figure = figure('color','w');
h{2}.axes = gca;
hgroup{2} = copyobj(hgroup{1},h{2}.axes);

% Formatting:
hold on
axis image
xlim([-1,1]*Lx_outer*0.7)
ylim([-1,1]*Ly_outer*0.7)
set(h{2}.axes,'YTick',[],'XTick',[],'XColor','w','YColor','w')

% Add white rectangle to trim image:
trimRectangle.x = [-1,0 ,0 ,-1,-1]*Lx;
trimRectangle.y = [-1,-1,+1,+1,-1]*Ly;
trimRectangle.color = 'w';

% "z" coordinate:
fields.String = '' ;
fields.Color = 'k';
x =      [0,1.5/100]*scaleFactor.x;
y =      [0,0    ];
h{2}.coordinate(1) = myTextArrow(gca,x,y,fields);
h{2}.coordinateText(1) = text(0.015,1.3e-4,'$\hat{z}_*$');

% "y" coordinate:
fields.String = '' ;
y =      [0,1.5/100]*scaleFactor.y;
x =      [0,0    ];
hta = myTextArrow(gca,x,y,fields);
h{2}.coordinateText(2) = text(0.1/100,2/100,'$\hat{y}_*$');

% Formatting coordinate:
set(h{2}.coordinateText,'interpreter','latex','FontSize',15)

% Trim rectangle:
htR = fill(trimRectangle.x,trimRectangle.y,trimRectangle.color);
set(htR,'EdgeColor',trimRectangle.color)

% Heat flux patter:
hG = plot(-heatFlux.sideView.z*Lz,heatFlux.sideView.y,'r');
set(hG,'lineWidth',2)
hText  = text(-0.0515, 0.0256,'$q(x_*,y_*,t)$');
set(hText,'interpreter','latex','FontSize',14,'color','r')

% Heat flux arrows:
fields.String = '' ;
fields.Color = 'r';
fields.LineWidth = 3;
for ss = 25:10:75
    y =      [1,1]*heatFlux.sideView.y(ss);
    x =      [-heatFlux.sideView.z(ss)*Lz*0.8,0];
    hta = myTextArrow(gca,x,y,fields);
end

% center of coordinate system:
plot(0,0,'ko','MarkerFaceColor','k')
 
% Dimensions 1:
% -------------------------------------------------------------------------
kk = 1;
fields.LineWidth = 2;
% thin line 1
nn = 1;
dim{kk}.thinline(nn).y1 = min(plate.y) - 0.2/100;
dim{kk}.thinline(nn).y2 = dim{kk}.thinline(nn).y1 - 1.5/100;
dim{kk}.thinline(nn).x1 = 0;
dim{kk}.thinline(nn).x2 = 0;
x1 = dim{kk}.thinline(nn).x1;
x2 = dim{kk}.thinline(nn).x2;
y1 = dim{kk}.thinline(nn).y1;
y2 = dim{kk}.thinline(nn).y2;
h{2}.dim{kk}.thinline(nn) = line([x1,x2],[y1,y2],'color','k');

% thin line 2
nn = 2;
dim{kk}.thinline(nn).y1 = min(plate.y) - 0.2/100;
dim{kk}.thinline(nn).y2 = dim{kk}.thinline(nn).y1 - 1.5/100;
dim{kk}.thinline(nn).x1 = max(plate.x);
dim{kk}.thinline(nn).x2 = max(plate.x);
x1 = dim{kk}.thinline(nn).x1;
x2 = dim{kk}.thinline(nn).x2;
y1 = dim{kk}.thinline(nn).y1;
y2 = dim{kk}.thinline(nn).y2;
h{2}.dim{kk}.thinline(nn) = line([x1,x2],[y1,y2],'color','k');

% Double arrow:
dim{kk}.doubleArrow.y1 = dim{kk}.thinline(nn).y2 + 0.2/100;
dim{kk}.doubleArrow.y2 = dim{kk}.doubleArrow.y1;
dim{kk}.doubleArrow.x1 = max(plate.x);
dim{kk}.doubleArrow.x2 = 0;

x1 = dim{kk}.doubleArrow.x1;
x2 = dim{kk}.doubleArrow.x2;
y1 = dim{kk}.doubleArrow.y1;
y2 = dim{kk}.doubleArrow.y2;
[figx,figy] = dsxy2figxy(gca,[x1 x2]*scaleFactor.x,[y1 y2]*scaleFactor.y);
h{2}.dim{kk}.darrow = annotation('doublearrow',figx,figy,'Color','k');
h{2}.dim{kk}.darrow.HeadStyle = 'vback2';

h{2}.dim{kk}.text{kk} = text(2/100,-0.073,'$L_z$');
set(h{2}.dim{kk}.text{kk},'interpreter','latex','FontSize',15)

% Save figure:
% =========================================================================
if saveFig
    figureName = 'Step_1_SideViewPlate';
    % PDF figure:
    exportgraphics(gcf,[figureName,'.pdf'],'Resolution',600); %,'ContentType', 'vector') 'BackgroundColor','none','ContentType','vector'

    % TIFF figure:
    exportgraphics(gcf,[figureName,'.tiff'],'Resolution',600) 
end

%% functions:
function circle = DrawCircle(p0,rad,N)
% p0: vector of center of circle
% r : radius of circle
% N: number of points

t = linspace(0,2*pi,N);
r.x = rad*cos(t);
r.y = rad*sin(t);

circle.x = p0.x + r.x;
circle.y = p0.y + r.y;
end