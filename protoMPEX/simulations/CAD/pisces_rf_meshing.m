close all
clear all

%Start Routine
nP = 1e6;

model = createpde;
importGeometry(model,'pisces_rf2.stl');% Import STL file
% importGeometry(model,'pisces-raf-geom_trunc4.stl');% Import STL file
figure(2)
pdegplot(model,'FaceLabels','on') %Plot stl 

tic %time meshing - for high resolution this is a large cost
mesh = generateMesh(model,'GeometricOrder','linear','Hmax',.1);% Options Hmax and Hmin can be set, linear order can also be used
%Hmax default was 0.005
figure(1)
pdeplot3D(model,'FaceAlpha',0.5)

[p,e,t] = meshToPet(mesh);

nPoints = length(p);
nTets = length(t);

tess = transpose(t(1:4,:));%sort(t(1:4,:),1);
% all faces
faces=[tess(:,[1 2 3]);tess(:,[1 2 4]); ...
       tess(:,[1 3 4]);tess(:,[2 3 4])];


faces = sort(faces,2);
faces = sortrows(faces);
Y = diff(faces);
zeroRow = [0,0,0];
k = ismember(Y,zeroRow,'rows');
k2 = find(k~=0);

faces([k2;k2+1],:) = [];

C = faces;

planes = zeros(length(C), 9);

planes(1:length(C),:) = [transpose(p(1:3,C(:,1))),transpose(p(1:3,C(:,2))),transpose(p(1:3,C(:,3)))];
toc

% tol = 1e-6;
% verticalInds1 = find(abs(planes(:,1) - planes(:,4)) < tol & abs(planes(:,1) - planes(:,7)) < tol);
% verticalInds2 = find(abs(planes(:,2) - planes(:,5)) < tol & abs(planes(:,2) - planes(:,8)) < tol);
% horizontalInds = find(abs(planes(:,3) - planes(:,6)) < tol & abs(planes(:,3) - planes(:,9)) < tol);
% walls = union(verticalInds1,verticalInds2);
% walls = union(walls,horizontalInds);
% allInds = 1:1:length(faces);
% allInds(walls) =[];
% materialSurfs = allInds;
% materialZ(allInds) = 74.0;


% plotSet = walls;
plotSet = 1:1:length(planes);
X = [planes((plotSet),1),planes((plotSet),4),planes((plotSet),7)];
Y = [planes((plotSet),2),planes((plotSet),5),planes((plotSet),8)];
Z = [planes((plotSet),3),planes((plotSet),6),planes((plotSet),9)];

% [X,Y,Z] = refineXYZ(X,Y,Z,6)
planes = [X(:,1) Y(:,1) Z(:,1) X(:,2) Y(:,2) Z(:,2) X(:,3) Y(:,3) Z(:,3)];
figure(10)
patch(transpose(X),transpose(Y),transpose(Z),'g','FaceAlpha',.3,'EdgeColor','k')%'none')
materialZ = 74*ones(length(planes),1);
surfs = ones(length(planes),1);

title({'PIESCES-RF Simulated GITR Geometry'})
xlabel('X [m]')
ylabel('Y [m]')
zlabel('Z [m]')
legend('Al')