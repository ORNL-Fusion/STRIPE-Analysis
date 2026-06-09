% close all
% clear all
% 
A = readmatrix('sheath.dat');
B = readmatrix('elements.dat')+1;

% AA=readtable('vsheath.dat');


i1 = A(B(:,2),3)';
i2 = A(B(:,3),3)';
i3 = A(B(:,4),3)';

j1 = A(B(:,2),4)';
j2 = A(B(:,3),4)';
j3 = A(B(:,4),4)';

k1 = A(B(:,2),5)';
k2 = A(B(:,3),5)';
k3 = A(B(:,4),5)';

p1 = A(B(:,2),6);
p2 = A(B(:,3),6);
p3 = A(B(:,4),6);

subset = 1:length(i1);

LL = [transpose(i1(subset)),transpose(i2(subset)),transpose(i3(subset))];
MM = [transpose(j1(subset)),transpose(j2(subset)),transpose(j3(subset))];
NN = [transpose(k1(subset)),transpose(k2(subset)),transpose(k3(subset))];

MM=-MM;
% p=sqrt(p1.^2+p2.^2+p3.^2)

figure

patch(transpose(LL),transpose(MM),transpose(NN),transpose(p1),'FaceAlpha',1,'EdgeAlpha', 0.3)
 hold on;


%%
if (exist('x1') == 0)

fid = fopen('gitrGeometryPointPlane3d.cfg');

tline = fgetl(fid);

tline = fgetl(fid);

    for i=1:18

        tline = fgetl(fid);

        evalc(tline);

    end

   Zsurface = Z;

end


hold on;



subset = 1:length(x1);%find(r<0.07 & z1> 0.001 & z1 < .20);


X = [transpose(x1(subset)),transpose(x2(subset)),transpose(x3(subset))];

Y = [transpose(y1(subset)),transpose(y2(subset)),transpose(y3(subset))];

Z = [transpose(z1(subset)),transpose(z2(subset)),transpose(z3(subset))];

%patch(transpose(X(surface,:)),transpose(Y(surface,:)),transpose(Z(surface,:)),impacts(surface),'FaceAlpha',.3)

% patch(transpose(X),transpose(Y),transpose(Z),zeros(1,length(subset)),'FaceAlpha',.3,'EdgeAlpha', 0.3)%,impacts(surface);


% theta_rotate = -0.167*pi;Xn = X.*cos(theta_rotate) - Y.*sin(theta_rotate);
% Yn = X.*sin(theta_rotate) + Y.*cos(theta_rotate);
% Zn = Z;
 
patch(transpose(X),transpose(Y),transpose(Z),zeros(1,length(subset)),'FaceAlpha',1,'EdgeAlpha', 1);hold on;
% patch(transpose(X),transpose(Y),transpose(Z),v_data,'FaceAlpha',1,'EdgeAlpha', 1)%,impacts(surface);
axis equal
