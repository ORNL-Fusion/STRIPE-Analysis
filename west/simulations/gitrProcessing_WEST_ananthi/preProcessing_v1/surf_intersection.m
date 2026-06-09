clc;close all; clear all;

% Read in GITR geometry file
% fid = fopen('/Users/78k/Library/CloudStorage/OneDrive-OakRidgeNationalLaboratory/ORNL-ATUL-MBP/myRepos/west/gitrProcessing_WEST/preProcessing_v1/gitrGeometryPointPlane3d_comsol.cfg');
fid = fopen('gitrGeometryPointPlane3d.cfg');
centroid=readmatrix('centroid.csv');
r_centroid=sqrt(centroid(:,1).^2+centroid(:,2).^2);

tline = fgetl(fid);
tline = fgetl(fid);
for i=1:18
    tline = fgetl(fid);
    evalc(tline);
end
subset = 1:length(x1);
X = [transpose(x1(subset)),transpose(x2(subset)),transpose(x3(subset))];
Y = [transpose(y1(subset)),transpose(y2(subset)),transpose(y3(subset))];
Z = [transpose(z1(subset)),transpose(z2(subset)),transpose(z3(subset))];
figure
patch(transpose(X),transpose(Y),transpose(Z),'g','FaceAlpha',0.4,'EdgeAlpha', 0.3)%,impacts(surface)
title('ProtoMPEX Geometry')
xlabel('x [m]')
ylabel('y [m]')
zlabel('z [m]')
set(gca,'fontsize',13)



p0 = [0 0 0];
% p1 = [2.99 2.125 0];
p1 = [0.36 3.9 -0.4]; 

% z_vertL = [0.4; 0.3473; 0.1689; 0.0749; 0.0106; -0.0851; -0.3042; -0.4];
% z_vertR = [0.4; 0.3713; 0.1828; 0.0839; 0.0034; -0.0913; -0.2782; -0.4];


intersected_surfaces = [];
ps=[];
for i=1:length(x1)
    % y1(i)=x1(i);
    % x1(i)=y1(i);
    % x1(i)=y1(i);

   [did_hit p]= surf_intersect(p0,p1,x1(i),y1(i),z1(i),x2(i),y2(i),z2(i),x3(i),y3(i),z3(i),a(i),b(i),c(i),d(i),plane_norm(i));
    if (did_hit)
        intersected_surfaces = [intersected_surfaces;i];
        ps=[ps;p];
        
    end
end

subset = intersected_surfaces;
X = [transpose(x1(subset)),transpose(x2(subset)),transpose(x3(subset))];
Y = [transpose(y1(subset)),transpose(y2(subset)),transpose(y3(subset))];
Z = [transpose(z1(subset)),transpose(z2(subset)),transpose(z3(subset))];
hold on;
patch(transpose(X),transpose(Y),transpose(Z),'b','FaceAlpha',1,'EdgeAlpha', 0.3)%,impacts(surface)
hold on;
plot3([p0(1) p1(1)],[p0(2) p1(2)],[p0(3) p1(3)],'LineWidth',2,'Color','r' )

r_centroid1=sqrt(centroid(intersected_surfaces,1).^2+centroid(intersected_surfaces,2).^2);
distance=vecnorm(ps')

% [value, index]=min(r_centroid1);
[value, index]=min(distance);
intersected_surfaces(index)
% ps(index)

% sxb=load('sxb_interp.csv');
% r=load('r.csv');
% z=('z.csv');
% 
% [m,n]=size(sxb);
% 
% sxb_new=interp2(sxb,3.0280/(max(r)-min(r))*m, 0.4/(max(z)-min(z))*n)
% 

function [did_hit p] = surf_intersect(p0,p1,x1,y1,z1,x2,y2,z2,x3,y3,z3,a,b,c,d,plane_norm)
did_hit = 0;
pointToPlaneDistance0 = (a * p0(1) + b * p0(2) + c * p0(3) + d) / plane_norm;
pointToPlaneDistance1 = (a * p1(1) + b * p1(2) + c * p1(3) + d) / plane_norm;

signPoint0 = sign(pointToPlaneDistance0);
signPoint1 = sign(pointToPlaneDistance1);

p=[NaN NaN NaN];

if (signPoint0 ~= signPoint1) 
      
        t = -(a * p0(1) + b * p0(2) + c * p0(3) + d) / ...
            (a * (p1(1) - p0(1)) + b * (p1(2) - p0(2)) + c * (p1(3) - p0(3)));

        p = [p0(1) + t * (p1(1) - p0(1)), p0(2) + t * (p1(2) - p0(2)), ...
                     p0(3) + t * (p1(3) - p0(3))];

A = [x1, y1,z1];
B = [x2, y2, z2];
C = [x3, y3, z3];
AB = B-A;
AC = C-A;
BC = C-B;
CA = A-C;

Ap = p-A;
Bp = p-B;
Cp = p-C;

normalVector = cross(AB,AC);
crossABAp = cross(AB,Ap);
crossBCBp = cross(BC,Bp);
crossCACp = cross(CA,Cp);


        signDot0 = sign(dot(crossABAp, normalVector));
        signDot1 = sign(dot(crossBCBp, normalVector));
        signDot2 = sign(dot(crossCACp, normalVector));
        totalSigns = abs(signDot0 + signDot1 + signDot2);



        hitSurface = 0;
        if (totalSigns == 3.0) 

          did_hit = 1;
          
        end
end
end