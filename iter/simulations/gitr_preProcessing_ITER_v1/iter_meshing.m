close all
clear all
[v, f, n, name] = stlRead('iter_extrusion.stl');
C = f;
p = v;

% planes is a variable which contains the xyz points of each triangle
% therefore it is a # triangles by 9 array
planes = zeros(length(C), 9);

planes(1:length(C),:) = [transpose(p(C(:,1),1:3))',transpose(p(C(:,2),1:3))',transpose(p(C(:,3),1:3))'];

plotSet= 1:length(C);
X = [planes((plotSet),1),planes((plotSet),4),planes((plotSet),7)];
Y = [planes((plotSet),2),planes((plotSet),5),planes((plotSet),8)];
Z = [planes((plotSet),3),planes((plotSet),6),planes((plotSet),9)];
    centroid = [1/3*(planes(:,1)+planes(:,4)+planes(:,7)), ...
        1/3*(planes(:,2)+planes(:,5)+planes(:,8)), ...
        1/3*(planes(:,3)+planes(:,6)+planes(:,9))];

%% Select regions for mesh refinement

disp('>>>> Refining geometry mesh')

limiter_inds = find((centroid(:,1) < 8500) & (abs(centroid(:,2))>800) & centroid(:,3) > -438 & centroid(:,3) < 1678);
% limiter_inds = find((abs(centroid(:,1)) < 1.5) & (centroid(:,2)>1.6) & centroid(:,3) > 0.2 & centroid(:,3) < 0.9);

plotSet= limiter_inds;
X = [planes((plotSet),1),planes((plotSet),4),planes((plotSet),7)];
Y = [planes((plotSet),2),planes((plotSet),5),planes((plotSet),8)];
Z = [planes((plotSet),3),planes((plotSet),6),planes((plotSet),9)];
figure(101)
patch(transpose(X),transpose(Y),transpose(Z),'g','FaceAlpha',1,'EdgeColor','k')

planes = [X(:,1) Y(:,1) Z(:,1) X(:,2) Y(:,2) Z(:,2) X(:,3) Y(:,3) Z(:,3)];
plotSet= 1:length(planes);
X = [planes((plotSet),1),planes((plotSet),4),planes((plotSet),7)];
Y = [planes((plotSet),2),planes((plotSet),5),planes((plotSet),8)];
Z = [planes((plotSet),3),planes((plotSet),6),planes((plotSet),9)];

centroid = [1/3*(planes(plotSet,1)+planes(plotSet,4)+planes(plotSet,7)), ...
        1/3*(planes(plotSet,2)+planes(plotSet,5)+planes(plotSet,8)), ...
        1/3*(planes(plotSet,3)+planes(plotSet,6)+planes(plotSet,9))];
%% Functions
function samples = sample_triangle(x,y,z,nP)
x_transform = x - x(1);
y_transform = y - y(1);
z_transform = z- z(1);

% figure(2)
% plot3([x_transform x_transform(1)],[y_transform y_transform(1)],[z_transform z_transform(1)])

v1 = [x_transform(2) y_transform(2) z_transform(2)];
v2 = [x_transform(3) y_transform(3) z_transform(3)];
v12 = v2 - v1;
normalVec = cross(v1,v2);

a1 = rand(nP,1);
a2 = rand(nP,1);

samples = a1.*v1 + a2.*v2;
% hold on
% scatter3(samples(:,1),samples(:,2),samples(:,3))
samples2x = samples(:,1) - v2(1);
samples2y = samples(:,2) - v2(2);
samples2z = samples(:,3) - v2(3);
samples12x = samples(:,1) - v1(1);
samples12y = samples(:,2) - v1(2);
samples12z = samples(:,3) - v1(3);
v1Cross = [(v1(2).*samples(:,3) - v1(3).*samples(:,2)) (v1(3).*samples(:,1) - v1(1).*samples(:,3)) (v1(1).*samples(:,2) - v1(2).*samples(:,1))];
v2 = -v2;
v2Cross = [(v2(2).*samples2z - v2(3).*samples2y) (v2(3).*samples2x - v2(1).*samples2z) (v2(1).*samples2y - v2(2).*samples2x)];
v12Cross = [(v12(2).*samples12z - v12(3).*samples12y) (v12(3).*samples12x - v12(1).*samples12z) (v12(1).*samples12y - v12(2).*samples12x)];

v1CD = normalVec(1)*v1Cross(:,1) + normalVec(2)*v1Cross(:,2) + normalVec(3)*v1Cross(:,3);
v2CD = normalVec(1)*v2Cross(:,1) + normalVec(2)*v2Cross(:,2) + normalVec(3)*v2Cross(:,3);
v12CD = normalVec(1)*v12Cross(:,1) + normalVec(2)*v12Cross(:,2) + normalVec(3)*v12Cross(:,3);

inside = abs(sign(v1CD) + sign(v2CD) + sign(v12CD));
insideInd = find(inside ==3);
notInsideInd = find(inside ~=3);
% scatter3(samples(insideInd,1),samples(insideInd,2),samples(insideInd,3))

 v2 = -v2;
dAlongV1 = v1(1).*samples(notInsideInd,1) + v1(2).*samples(notInsideInd,2) + v1(3).*samples(notInsideInd,3);
dAlongV2 = v2(1).*samples(notInsideInd,1) + v2(2).*samples(notInsideInd,2) + v2(3).*samples(notInsideInd,3);

dV1 = norm(v1);
dV2 = norm(v2);
halfdV1 = 0.5*dV1;
halfdV2 = 0.5*dV2;

samples(notInsideInd,:) = [-(samples(notInsideInd,1) - 0.5*v1(1))+0.5*v1(1) ...
    -(samples(notInsideInd,2) - 0.5*v1(2))+0.5*v1(2) ...
    -(samples(notInsideInd,3) - 0.5*v1(3))+0.5*v1(3)];
% samples(notInsideInd,:) = [-(samples(notInsideInd,1) - 0.5*v2(1))+0.5*v2(1) ...
%     -(samples(notInsideInd,2) - 0.5*v2(2))+0.5*v2(2) ...
%     -(samples(notInsideInd,3) - 0.5*v2(3))+0.5*v2(3)];
samples(notInsideInd,:) = [(samples(notInsideInd,1) + v2(1)) ...
    (samples(notInsideInd,2) +v2(2)) ...
    (samples(notInsideInd,3) + v2(3))];
% figure(4)
% plot3([x_transform x_transform(1)],[y_transform y_transform(1)],[z_transform z_transform(1)])
% hold on
% scatter3(samples(:,1),samples(:,2),samples(:,3))

samples(:,1) = samples(:,1)+ x(1);
samples(:,2) = samples(:,2)+ y(1);
samples(:,3) = samples(:,3)+ z(1);

% figure(5)
% plot3([x x(1)],[y y(1)],[z z(1)])
% hold on
% scatter3(samples(:,1),samples(:,2),samples(:,3))
end




function [planes, tet_centers] = get_triangles_from_stl(file_string,color,resolution,scale)
model = createpde;
importGeometry(model,file_string);
figure(1)
pdegplot(model,'FaceLabels','on') %Plot stl

tic %time meshing - for high resolution this is a large cost
mesh = generateMesh(model,'GeometricOrder','linear','Hmax',resolution);% Options Hmax and Hmin can be set, linear order can also be used

figure(2)
pdeplot3D(model,'FaceAlpha',0.5)

% Convert mesh to PET format (points, edges, triangles)
[p,e,t] = meshToPet(mesh);
p_copy=p;
p(2,:)=p_copy(3,:);
p(3,:)=p_copy(2,:);

nPoints = length(p);
nTets = length(t);

% Take 3D tetrahedral mesh and perform "skinning" operation
% Skinning operation takes the "outward facing triangles" from the
% mesh. This is used as the 3D surface mesh for GITR
tess = transpose(t(1:4,:));

% all faces
faces=[tess(:,[1 2 3]);tess(:,[1 2 4]); ...
    tess(:,[1 3 4]);tess(:,[2 3 4])];

% inds = ceil(0:0.25:(length(tess)-0.25));



inds = [1:length(tess) 1:length(tess) 1:length(tess) 1:length(tess)];


faces = sort(faces,2);
[faces, i] = sortrows(faces);



inds = inds(i);

Y = diff(faces);
zeroRow = [0,0,0];
k = ismember(Y,zeroRow,'rows');
k2 = find(k~=0);

faces([k2;k2+1],:) = [];
inds([k2;k2+1]) = [];
tet_center_x = 0.25*sum([p(1,tess(inds,1))',p(1,tess(inds,2))',p(1,tess(inds,3))',p(1,tess(inds,4))'],2);
tet_center_y = 0.25*sum([p(2,tess(inds,1))',p(2,tess(inds,2))',p(2,tess(inds,3))',p(2,tess(inds,4))'],2);
tet_center_z = 0.25*sum([p(3,tess(inds,1))',p(3,tess(inds,2))',p(3,tess(inds,3))',p(3,tess(inds,4))'],2);
tet_centers = [tet_center_x,tet_center_y,tet_center_z];
% end of skinning operation

C = faces;

% planes is a variable which contains the xyz points of each triangle
% therefore it is a # triangles by 9 array
planes = zeros(length(C), 9);

planes(1:length(C),:) = [transpose(p(1:3,C(:,1))),transpose(p(1:3,C(:,2))),transpose(p(1:3,C(:,3)))]./scale;
toc


plotSet = 1:1:length(planes);
X = [planes((plotSet),1),planes((plotSet),4),planes((plotSet),7)];
Y = [planes((plotSet),2),planes((plotSet),5),planes((plotSet),8)];
Z = [planes((plotSet),3),planes((plotSet),6),planes((plotSet),9)];

planes = [X(:,1) Y(:,1) Z(:,1) X(:,2) Y(:,2) Z(:,2) X(:,3) Y(:,3) Z(:,3)];
figure(3)
hold on
patch(transpose(X),transpose(Y),transpose(Z),color,'FaceAlpha',1,'EdgeColor','k')%'none')
end