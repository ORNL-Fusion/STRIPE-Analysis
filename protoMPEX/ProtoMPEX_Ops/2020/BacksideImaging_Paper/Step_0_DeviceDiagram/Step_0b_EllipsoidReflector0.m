% Step_0b:
% =========================================================================
% Ray tracing:

clear all
close all
clc

saveFig = 0;

% Conversion factor:
% =========================================================================
in2m = 2.54/100;

% Create reflection surface:
% ==========================
% Ellipse parameters:
surface.a  = 0.23467;
surface.b  = 0.18816;
surface.F0 = [-0.04664,0.24252]';
surface.F1 = [-0.0762 ,0      ]';

% Ellipse origin vector:
% ======================
surface.F3 = (surface.F0 + surface.F1)/2;

% Ellipse alliged coordinate system:
% =================================
% e1:
dF = surface.F0 - surface.F1;
e1 = dF/norm(dF);
% e2:
a1 = e1(1);
b1 = e1(2);
b2 = 1/sqrt( 1 + (b1/a1)^2);
a2 = -(b1/a1)*b2;
e2 = [a2,b2]';
% Coordinate system matrix:
e  = [e1,e2];

% Ellipse equation on "e" coordinate system:
% ==========================================
x_e = linspace(-surface.a,surface.a,5e3)';
y_e = -surface.b*sqrt(1 - (x_e/surface.a).^2);
p_e = [x_e,y_e]';

% Transform to "standard" coordinate system:
% =========================================================================
% Standard coordinate system:
s1 = [1 0]';
s2 = [0 1]';
s  = [s1,s2];

% Transform to "s" coordinate system:
p_s = e*p_e;

% Ellipse position vector on "s" coordinate system:
m_s = p_s + surface.F3;

% Select ellipse segment:
% =========================================================================
yMin = -2.5*2.54/100;
yMax = -yMin;
rng_y = find(m_s(2,:) >= yMin & m_s(2,:) <= yMax);

% Create surface vector:
% =========================================================================
surface.x = m_s(1,rng_y)';
surface.y = m_s(2,rng_y)';

% Add surface base:
% =========================================================================
x1 = surface.x(1);
x2 = 0 - 0*0.625*in2m;
x3 = x2 + 2.75*in2m;
x4 = x2 + 5.68*in2m;
x5 = surface.x(end);

y1 = surface.y(1);
y2 = y1 - 0.39*in2m;
y3 = surface.y(end);

base.x = [x5, x4, x3, x2, x2, x1]';
base.y = [y3, y3, y2, y2, y1, y1]';

% Assemble mirror:
% =========================================================================
mirror.x = [ base.x; surface.x];
mirror.y = [ base.y; surface.y];

% Offset relative to chamber:
% =========================================================================
rho.x = 0;
rho.y = -12.315*in2m/2;

% Create plasma cross section:
% =========================================================================
p0.x = 0;
p0.y = 0;
r = 0.04;
plasma = DrawCircle(p0,r,500);

% p0.x = -0.6;
% p0.y = 0.6;
% r = 0.7;
% plasma = DrawCircle(p0,r,500);

% Target plate boundary:
% =========================================================================
target.width = 10/100;
target.height = 10/100;
x1 = target.width/2;
y1 = target.width/2;
t_x = [+x1,-x1,-x1,+x1,+x1];
t_y = [+y1,+y1,-y1,-y1,+y1];
R = [cosd(45) sind(45); -sind(45) cosd(45)];
tt = R*[t_x;t_y];
target.x = tt(1,:);
target.y = tt(2,:);

% Create chamber cross section:
% =========================================================================
p0.x = 0;
p0.y = 0;
r = 19.885*in2m/2;
vacBoundary1 = DrawCircle(p0,r,500);

p0.x = 0;
p0.y = 0;
r = 22.25*in2m/2;
vacBoundary2 = DrawCircle(p0,r,500);

% Create waveguide:
% =========================================================================
waveguide.ID = 63.5/1000;
x1 = -0.3 ;
x2 = surface.F1(1) ;
y1 = surface.F1(2) + waveguide.ID/2;

waveguide.x = [x1, x2 , x2, x1];
waveguide.y = [y1, y1 ,-y1,-y1];

% Create reflection surface:
% =========================================================================
reflectiveSurface.x = surface.x + rho.x;
reflectiveSurface.y = surface.y + rho.y;

% Create absorption surface:
% =========================================================================
rng = find(vacBoundary1.y > 0 );
absorptionSurface.x = vacBoundary1.x(rng);
absorptionSurface.y = vacBoundary1.y(rng);

rng = find(plasma.x > 0 & plasma.y < 0);
absorptionSurface.x = plasma.x(rng);
absorptionSurface.y = plasma.y(rng);

% rng = find(plasma.x < 0);
% absorptionSurface.x = plasma.x(rng);
% absorptionSurface.y = plasma.y(rng);

% Create ray:
% =========================================================================
n_ray = 15;
dy     = linspace(-0.9*waveguide.ID/2,0.9*waveguide.ID/2,n_ray);
dtheta = linspace(-7,7,n_ray);

for jj = 1:numel(dy)
    ray{jj}.x(1)  = waveguide.x(2) + rho.x;
    ray{jj}.y(1)  = surface.F1(2) + rho.y + dy(jj);
    ray{jj}.theta(1) = dtheta(jj);
    
    ray{jj} = intersectPoint(ray{jj},reflectiveSurface);
    ray{jj} = intersectPoint(ray{jj},absorptionSurface);
end

%% Plot ellipse segment:
% =========================================================================
close all

figure; 
hold on
plot(surface.F0(1),surface.F0(2),'k.')
plot(surface.F1(1),surface.F1(2),'k.')
plot(surface.F3(1),surface.F3(2),'k.')
plot(mirror.x,mirror.y,'k-')
plot(mirror.x,mirror.y,'k-')
plot(m_s(1,:),m_s(2,:),'r-')

axis image
xlim([-0.3,0.3])
ylim([-0.3,0.4])
grid on

figure('color','w'); 
showFoci = 0;
showTarget = 1;
box on
hold on
if showFoci
    plot(surface.F0(1)+ rho.x,surface.F0(2) + rho.y,'go')
    plot(surface.F1(1)+ rho.x,surface.F1(2) + rho.y,'go')
end
if showTarget
    plot(target.x ,target.y,'k--','lineWidth',1)
end
plot(mirror.x + rho.x,mirror.y + rho.y,'k-','lineWidth',2)
plot(plasma.x,plasma.y,'k.','lineWidth',2)
plot(vacBoundary1.x,vacBoundary1.y,'k-','lineWidth',2)
% plot(vacBoundary2.x,vacBoundary2.y,'k-','lineWidth',2)
plot(waveguide.x + rho.x ,waveguide.y + rho.y, 'k-','lineWidth',2)
for jj = 1:numel(dy)
    plot(ray{jj}.x,ray{jj}.y,'g','lineWidth',1)
end

% Formatting:
axis square
xlim([-0.3,0.3])
ylim([-0.3,0.3])
set(gca,'FontName','Times','FontSize',11)
xlabel('x [m]','interpreter','Latex','FontSize',13)
ylabel('y [m]','interpreter','Latex','FontSize',13)
title(['Sectional-view A-A'],'Interpreter','Latex','FontSize',13)

% Arrows:
% =========================================================================
% Text arrows fields:
fields.HorizontalAlignment = 'center';
fields.Color = 'k';
fields.FontSize = 11;
fields.Interpreter = 'Latex';
fields.HeadLength = 5;
fields.HeadWidth = 5;
fields.HeadStyle = 'vback2';

% Arrow convention:
% [end,start]

% 28 GHz waveguide:
fields.String = 'D';
x =      [-0.12 ,-0.12];
y =      [-0.08 ,-0.125];
hta = myTextArrow(gca,x,y,fields);

% vacuum boundary:
fields.String = 'E';
x =      [-0.12 ,-0.18];
y =      [+0.1 ,+0.1];
hta = myTextArrow(gca,x,y,fields);

% Target boundary:
fields.String = 'F';
x =      [-0   ,-0    ];
y =      [+0.1 ,+0.07];
hta = myTextArrow(gca,x,y,fields);

% Plasma boundary:
fields.String = 'B';
x =      [+0.05,+0.025];
y =      [+0.05,+0.025];
hta = myTextArrow(gca,x,y,fields);

% Ellipsoid reflector:
fields.String = 'G';
x =      [+0.1 ,+0.1];
y =      [-0.03,-0.09];
hta = myTextArrow(gca,x,y,fields);

% Save figure:
% =========================================================================
if saveFig
    figureName = 'Step_0b_SectionAA';
% PDF figure:
exportgraphics(gcf,[figureName,'.pdf'],'Resolution',600,'ContentType', 'vector') 

% TIFF figure:
exportgraphics(gcf,[figureName,'.tiff'],'Resolution',600) 
end

%% Functions:
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

function ray = intersectPoint(ray,surface)
% Ray conditions:
n  = numel(ray.y);
m0 = tand(ray.theta(n));

% Starting location of input ray:
y0 = ray.y(n);
x0 = ray.x(n);

% Components of the unit vector of input ray:
cos0 = 1/sqrt(1 + m0^2);
sin0 = m0/sqrt(1 + m0^2);

% Input unit vector:
r0 = [cos0;sin0];

% Surface conditions:
mi = diff(surface.y)./diff(surface.x);
xi = surface.x;
yi = surface.y;

% Calculate intersection points:
for ii = 1:(numel(mi) - 1)
    MI = [-m0 1; -mi(ii) 1];
    r  = inv(MI)*[y0 - m0*x0; yi(ii) - mi(ii)*xi(ii)];
    xI(ii) = r(1);
    yI(ii) = r(2);    
end

for ii = 1:(numel(mi) - 1)
    if ( xI(ii) >= xi(ii) && xI(ii) < xi(ii+1) )
       
        ray.x(n+1) = xI(ii);
        ray.y(n+1) = yI(ii);
        
        % Slope of intersecting segment:
        m1 = mi(ii);
        
        % Components of the unit vector along intersecting segment:
        cos1 = 1/sqrt(1 + m1^2);
        sin1 = m1/sqrt(1 + m1^2);
      
        % Rotation matrix:
        R1 = [+cos1, -sin1 ;...
              +sin1, +cos1];
        
        % Express ray unit vector in terms of segment's coordinate system:
        r1 = transpose(R1)*r0;
        
        % Reflect "y'" component of unit vector:
        r1(2) = -r1(2);
        
        % Transform back to standard coordinate system:
        r2 = R1*r1;
        
        % Calculate slope of reflecting ray:
        ray.theta(n+1) = atan2d(r2(2),r2(1));
        break
    end
end
end
