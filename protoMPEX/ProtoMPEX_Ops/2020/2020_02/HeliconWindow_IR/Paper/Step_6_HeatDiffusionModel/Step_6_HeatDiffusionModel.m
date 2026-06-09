% Step 6:
% Create the 1D heat diffusion model schematic:

clear all
clc
close all

figureName = 'Step_6_HeatDiffusionModel';

% Load slab image:
% =========================================================================
im = imread('Slab_1D.jpg');

% Create coordinates for the image:
% =========================================================================
[Ny,Nx,~] = size(im);
Lx = 2;
Ly = 1;
yy = linspace(-Ly,+Ly)/2;
xx = linspace(+00,+Lx);

% Conditions for arrows representing heat flux vectors:
% =========================================================================
dy = Ly/10;
nArrow = 10;
xL = Lx/5;
xGap = Lx/20;
x0 = min(xx) -xL - xGap;
yP = linspace (min(yy),max(yy),nArrow);
xP = [x0,x0 + xL];

%%
% Plot image:
% =========================================================================
figure('color','w')
hold on

% Slab:
him = image(xx,yy,im);
% Set xlimit:
xlim([-1.5,2.5])
ylim([-1,1])

% Arrows to represent heat flux:
for ii = 1:nArrow
    [xArrow,yArrow] = dsxy2figxy([xP(1) xP(2)],[yP(ii) yP(ii)]);
    hArrow(ii) = annotation('arrow',xArrow,yArrow);
end
set(hArrow,'HeadStyle','vback1','LineWidth',2,'color','r')

% Draw coordinate system:
fontSize1 = 14;
xArrow = [0 1];
yArrow = [0 0];
[xdum,ydum] = dsxy2figxy(xArrow,yArrow);
hUnitVec = annotation('arrow',xdum,ydum);
set(hUnitVec,'HeadStyle','vback1','LineWidth',1,'color','k')
text(xArrow(2) + Lx/20 ,yArrow(1),'$x$','Interpreter','Latex','fontsize',fontSize1)

% Length arrow:
xArrow = [min(xx), max(xx) - Lx/8];
yArrow = max(yy)*[1 1] + Ly/10;
[xdum,ydum] = dsxy2figxy(xArrow,yArrow);
hLengthArrow = annotation('doublearrow',xdum,ydum);
set(hLengthArrow,'HeadStyle','vback1','LineWidth',0.1,'color','k')
hLengthArrow.LineWidth = 0.1;
hLengthArrow.Head1Length = 7;
hLengthArrow.Head2Length = 7;
hLengthArrow.Head1Width = 5;
hLengthArrow.Head2Width = 5;
text(mean(xx)-Lx/10,max(yy)+Ly/6,'$L$','Interpreter','Latex','fontsize',fontSize1)

% Annotations:
fontSize2 = 15;
fontSize3 = 16;
% Temperature:
text(mean(xx) - Lx/7, Ly/4,'$T(x,t)$','Interpreter','Latex','fontsize',fontSize2)
% Heat flux:
text(xP(1) - Lx/5,0,'$q(t)$','Interpreter','Latex','fontsize',fontSize2)
% BC1:
text(xP(1) - Lx/8,min(yP) - Ly/5,'$-k\frac{\partial T}{\partial x}|_{x=0}=q(t)$','Interpreter','Latex','fontsize',fontSize3)
% BC2:
text(xP(2) + 0.66*Lx,min(yP) - Ly/5,'$-k\frac{\partial T}{\partial x}|_{x=L}=0$','Interpreter','Latex','fontsize',fontSize3)

% Final formatting:
set(gca,'Visible','off')

% Saving figure:
% =========================================================================
saveas(gcf,figureName,'tiffn')