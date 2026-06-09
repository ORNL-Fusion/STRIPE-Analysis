% close all
% clear all

% return
read_efit_data;
% % ITER_antenna_data
antenna_ITER_test_new;


figure
patch(transpose(X),transpose(Y),transpose(Z),'k','FaceAlpha',1,'EdgeColor','k')
scatter3(centroid(:,1),centroid(:,2),centroid(:,3));

%% --- Load Data ---
comsol_data = readmatrix('../comsol/1540_i3_VDC4.txt');
xx = comsol_data(:,1);
yy = comsol_data(:,2);
zz = comsol_data(:,3);
emag = comsol_data(:,4); % Field magnitude

% %% --- Method 1: 3D Scatter Plot ---
% figure;
% scatter3(xx, yy, zz, 15, emag, 'filled'); % '15' controls marker size
% colormap parula;
% colorbar;
% title('3D Electric Field Magnitude');
% xlabel('X [m]');
% ylabel('Y [m]');
% zlabel('Z [m]');
% view(3); % Set 3D view
% grid on;
% set(gca, 'FontSize', 14, 'FontName', 'Times');
% 
% %% --- Method 2: Smooth 3D Surface Plot (If Data is Structured) ---
% tri = delaunayTriangulation(xx, yy, zz); % Triangulation for surface
% figure;
% trisurf(tri.ConnectivityList, xx, yy, zz, emag, 'EdgeColor', 'none');
% 
% % Apply smooth shading
% shading interp;
% lighting phong;
% material dull;
% 
% colormap parula;
% colorbar;
% title('3D Field Distribution (Smoothed)');
% xlabel('X [m]');
% ylabel('Y [m]');
% zlabel('Z [m]');
% view(3);
% axis equal;
% set(gca, 'FontSize', 14, 'FontName', 'Times');
figure; scatter3(xx,yy,zz,[], emag);

for i=1:length(centroid)
    distance = sqrt((centroid(i,1) - xx).^2 + (centroid(i,2) - yy).^2 + (centroid(i,3) - zz).^2);

    [M I] = min(distance);
    potential_surf(i) = emag(I);
end
potential_surf=potential_surf';


% % Write surface data
% 
surface_variables=table(potential_surf,ne_surf,te_surf,ti_surf,vNeMag_surf,Bmag_surf,Br_surf,Bt_surf,Bz_surf,theta,nNe_surf);
writetable(surface_variables,'Targets.txt')
figure
patch(transpose(X),transpose(Y),transpose(Z),'k','FaceAlpha',1,'EdgeColor','k');hold on;
scatter3(xx,yy,zz);

title('ITER Limiter Geometry For GITR')        
xlabel('x [mm]')
ylabel('y [mm]')
zlabel('z [mm]')

figure;
patch(transpose(X),transpose(Y),transpose(Z),ne_surf,'FaceAlpha',1,'EdgeAlpha', 0.3 );
title('ITER Limiter surface density')        
xlabel('x [mm]')
ylabel('y [mm]')
zlabel('z [mm]')

figure; 
patch(transpose(X),transpose(Y),transpose(Z),potential_surf,'FaceAlpha',1,'EdgeAlpha', 0.3 )
title('ITER Limiter surface potential')        
xlabel('x [mm]')
ylabel('y [mm]')
zlabel('z [mm]')



disp('>>>> Generatting gmtry file for GITR')

format = '%.16e';

fileID = fopen('gitrGeometryPointPlane3d.cfg','w');
fprintf(fileID,'geom = \n{ \n   x1 = [');
fprintf(fileID,format,planes(1,1));
for i=2:nFaces
    fprintf(fileID, ',');
    fprintf(fileID,format,planes(i,1));
end
fprintf(fileID,' ] \n   y1 = [');
fprintf(fileID,format,planes(1,2));
for i=2:nFaces
    fprintf(fileID, ',');
    fprintf(fileID,format,planes(i,2));
end
fprintf(fileID,' ] \n   z1 = [');
fprintf(fileID,format,planes(1,3));
for i=2:nFaces
    fprintf(fileID, ',');
    fprintf(fileID,format,planes(i,3));
end
fprintf(fileID,' ] \n   x2 = [');
fprintf(fileID,format,planes(1,4));
for i=2:nFaces
    fprintf(fileID, ',');
    fprintf(fileID,format,planes(i,4));
end
fprintf(fileID,' ] \n   y2 = [');
fprintf(fileID,format,planes(1,5));
for i=2:nFaces
    fprintf(fileID, ',');
    fprintf(fileID,format,planes(i,5));
end
fprintf(fileID,' ] \n   z2 = [');
fprintf(fileID,format,planes(1,6));
for i=2:nFaces
    fprintf(fileID, ',');
    fprintf(fileID,format,planes(i,6));
end
fprintf(fileID,' ] \n   x3 = [');
fprintf(fileID,format,planes(1,7));
for i=2:nFaces
    fprintf(fileID, ',');
    fprintf(fileID,format,planes(i,7));
end
fprintf(fileID,' ] \n   y3 = [');
fprintf(fileID,format,planes(1,8));
for i=2:nFaces
    fprintf(fileID, ',');
    fprintf(fileID,format,planes(i,8));
end
fprintf(fileID,' ] \n   z3 = [');
fprintf(fileID,format,planes(1,9));
for i=2:nFaces
    fprintf(fileID, ',');
    fprintf(fileID,format,planes(i,9));
end
fprintf(fileID,' ] \n   a = [');
fprintf(fileID,format,abcd(1,1));
for i=2:nFaces
    fprintf(fileID, ',');
    fprintf(fileID,format,abcd(i,1));
end

fprintf(fileID,' ] \n   b = [');
fprintf(fileID,format,abcd(1,2));
for i=2:nFaces
    fprintf(fileID, ',');
    fprintf(fileID,format,abcd(i,2));
end

fprintf(fileID,' ] \n   c = [');
fprintf(fileID,format,abcd(1,3));
for i=2:nFaces
    fprintf(fileID, ',');
    fprintf(fileID,format,abcd(i,3));
end
fprintf(fileID,' ] \n   d = [');
fprintf(fileID,format,abcd(1,4));
for i=2:nFaces
    fprintf(fileID, ',');
    fprintf(fileID,format,abcd(i,4));
end
fprintf(fileID,' ] \n   plane_norm = [');
fprintf(fileID,format,plane_norm(1));
for i=2:nFaces
    fprintf(fileID, ',');
    fprintf(fileID,format,plane_norm(i));
end
% fprintf(fileID,' ] \n   ABxAC = [');
% fprintf(fileID,format,ABxAC(1))
% for i=2:nFaces
% fprintf(fileID, ',')
% fprintf(fileID,format,ABxAC(i))
% end
fprintf(fileID,' ] \n   BCxBA = [');
fprintf(fileID,format,BCxBA(1));
for i=2:nFaces
    fprintf(fileID, ',');
    fprintf(fileID,format,BCxBA(i));
end
fprintf(fileID,' ] \n   CAxCB = [');
fprintf(fileID,format,CAxCB(1));
for i=2:nFaces
    fprintf(fileID, ',');
    fprintf(fileID,format,CAxCB(i));
end
fprintf(fileID,' ] \n   area = [');
fprintf(fileID,format,area(1,1));
for i=2:nFaces
    fprintf(fileID, ',');
    fprintf(fileID,format,area(i));
end
fprintf(fileID,' ] \n   Z = [');
fprintf(fileID,'%f',materialZ(1));
for i=2:nFaces
    fprintf(fileID, ',');
    fprintf(fileID,'%f',materialZ(i));
end
fprintf(fileID,' ] \n   surface = [');
fprintf(fileID,'%i',surfs(1));
for i=2:nFaces
    fprintf(fileID, ',');
    fprintf(fileID,'%i',surfs(i));
end
fprintf(fileID,' ] \n   inDir = [');
fprintf(fileID,'%i',inDir(1));
for i=2:nFaces
    fprintf(fileID, ',');
    fprintf(fileID,'%i',inDir(i));
end

fprintf(fileID,' ] \n   potential = [');
fprintf(fileID,'%f',potential_surf(1));
for i=2:nFaces
    fprintf(fileID, ',');
    fprintf(fileID,'%f',potential_surf(i));
end

fprintf(fileID,' ] \n');
fprintf(fileID,'periodic = 0;\n');
fprintf(fileID,'theta0 = 0.0;\n');
fprintf(fileID,'theta1 = 0.0\n');
fprintf(fileID,'periodic_bc_x0 = 0.0;\n');
fprintf(fileID,'periodic_bc_x1 = 0.0;\n');
fprintf(fileID,'periodic_bc_x = 0;}\n');


fclose(fileID);


%% Particle intialization
disp('>>>> Initializating Particles')

surf_inds = find(surfs);
nP = 1000;
nR = 100;

%% === Load Input Data ===
dataRF = readmatrix('Targets.txt');

% Load sputtering yields for Neon charge states
yields_Ne1 = readmatrix('yields_Ne+.csv');
yields_Ne2 = readmatrix('yields_Ne2+.csv');
yields_Ne3 = readmatrix('yields_Ne3+.csv');
yields_Ne4 = readmatrix('yields_Ne3+.csv');
yields_Ne5 = readmatrix('yields_Ne5+.csv');
yields_Ne6 = readmatrix('yields_Ne6+.csv');
yields_Ne7 = readmatrix('yields_Ne7+.csv');

% Neon fractional contributions (example values, update if needed)
ni_dataRF_Ne1 = 0.02 * 0.23 * dataRF(:,11);
ni_dataRF_Ne2 = 0.02 * 0.26 * dataRF(:,11);
ni_dataRF_Ne3 = 0.02 * 0.205 * dataRF(:,11);
ni_dataRF_Ne4 = 0.02 * 0.12 * dataRF(:,11);
ni_dataRF_Ne5 = 0.02 * 0.067 * dataRF(:,11);
ni_dataRF_Ne6 = 0.02 * 0.031 * dataRF(:,11);
ni_dataRF_Ne7 = 0.02 * 0.013 * dataRF(:,11);

% Velocity
v_dataRF_Ne = abs(dataRF(:,5));

potential_data=dataRF(:,1);
ne_data=dataRF(:,2);
te_data=dataRF(:,3);
v_data=dataRF(:,5);
erosion_flux=[0;yields_Ne5(:,5)].*ni_dataRF_Ne5.*v_dataRF_Ne;
erosion_flux(isnan(erosion_flux))=0;

% Erosion rate (Ne1+) = 1.0636e15 particles/s
% Erosion rate (Ne2+) = 1.2038e15 particles/s
% Erosion rate (Ne3+) = 9.4873e14 particles/s
% Erosion rate (Ne4+) = 5.5535e14 particles/s
% Erosion rate (Ne5+) = 3.0963e14 particles/s
% Erosion rate (Ne6+) = 1.4328e14 particles/s
% Erosion rate (Ne6+) = 1.4328e14 particles/s
% Erosion rate (Ne6+) = 0.60153e14 particles/s




% Erosion data
% Y0 = readmatrix('yields.csv');
% eroded_flux= ; % proportional to density now, need some realistic values
erosion = erosion_flux.*area;

erosion=erosion(surf_inds);
erosion_inds = find(erosion);
erosion_sub = erosion(erosion_inds);
erosion_sub_cdf = cumsum(erosion_sub);
erosion_rate = erosion_sub_cdf(end);
erosion_sub_cdf = erosion_sub_cdf./erosion_sub_cdf(end);


rand1 = rand(nP,1);

element = interp1([0, erosion_sub_cdf'],0:1:length(erosion_sub_cdf),rand1);

element_ceil = ceil(element);

% Particle initialization
xP = zeros(1,nP);
yP = zeros(1,nP);
zP = zeros(1,nP);
vxP = zeros(1,nP);
vyP = zeros(1,nP);
vzP = zeros(1,nP);


offset=1e-5;
for j=1:nP
    i = erosion_inds(element_ceil(j));
    normal = -abcd(i,1:3)/plane_norm(i);
    normal = inDir(i)*normal;

    x_tri = X(i,:)+offset*normal(1);
    y_tri = Y(i,:)+offset*normal(2);
    z_tri = Z(i,:)+offset*normal(3);

    samples = sample_triangle(x_tri,y_tri,z_tri,1);
    xP(j) = samples(:,1);
    yP(j) = samples(:,2);
    zP(j) = samples(:,3);
    vxP(j) = 5000*normal(1); %[m/s]
    vyP(j) = 5000*normal(2); %[m/s]
    vzP(j) = 5000*normal(3); %[m/s]
end
figure(301);
patch(transpose(X),transpose(Y),transpose(Z),'b','FaceAlpha',.3,'EdgeColor','none')%'none')
hold on
scatter3(xP,yP,zP,'r')

ncid = netcdf.create('./particle_source_iter.nc','NC_WRITE');

dimP = netcdf.defDim(ncid,'nP',nP);
xVar = netcdf.defVar(ncid,'x','double',[dimP]);
yVar = netcdf.defVar(ncid,'y','double',[dimP]);
zVar = netcdf.defVar(ncid,'z','double',[dimP]);
vxVar = netcdf.defVar(ncid,'vx','double',[dimP]);
vyVar = netcdf.defVar(ncid,'vy','double',[dimP]);
vzVar = netcdf.defVar(ncid,'vz','double',[dimP]);

netcdf.endDef(ncid);

netcdf.putVar(ncid, xVar, xP);
netcdf.putVar(ncid, yVar, yP);
netcdf.putVar(ncid, zVar, zP);
netcdf.putVar(ncid, vxVar, vxP);
netcdf.putVar(ncid, vyVar, vyP);
netcdf.putVar(ncid, vzVar, vzP);

netcdf.close(ncid);

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


function [Xrefined, Yrefined, Zrefined] = refineXYZ(X,Y,Z,n)

for j=1:n
    Xrefined = zeros(2*length(X),3);
    Yrefined = zeros(2*length(X),3);
    Zrefined = zeros(2*length(X),3);
    for i=1:length(X)
        A = [X(i,1) Y(i,1) Z(i,1)];
        B = [X(i,2) Y(i,2) Z(i,2)];
        C = [X(i,3) Y(i,3) Z(i,3)];

        AB = B-A;
        AC = C-A;
        BC = C-B;

        norms =[norm(AB) norm(BC) norm(AC)];
        [maxVal maxInd] = max(norms);
        if maxInd ==1
            midPtAB = A + 0.5*AB;
            Xrefined(2*i-1,1) = A(1);
            Xrefined(2*i,1) = midPtAB(1);
            Xrefined(2*i-1,2) = midPtAB(1);
            Xrefined(2*i,2) = B(1);
            Xrefined(2*i-1,3) = C(1);
            Xrefined(2*i,3) = C(1);
            Yrefined(2*i-1,1) = A(2);
            Yrefined(2*i,1) = midPtAB(2);
            Yrefined(2*i-1,2) = midPtAB(2);
            Yrefined(2*i,2) = B(2);
            Yrefined(2*i-1,3) = C(2);
            Yrefined(2*i,3) = C(2);
            Zrefined(2*i-1,1) = A(3);
            Zrefined(2*i,1) = midPtAB(3);
            Zrefined(2*i-1,2) = midPtAB(3);
            Zrefined(2*i,2) = B(3);
            Zrefined(2*i-1,3) = C(3);
            Zrefined(2*i,3) = C(3);
        elseif maxInd ==2
            midptBC = B + 0.5*BC;
            Xrefined(2*i-1,1) = A(1);
            Xrefined(2*i,1) = A(1);
            Xrefined(2*i-1,2) = B(1);
            Xrefined(2*i,2) = midptBC(1);
            Xrefined(2*i-1,3) = midptBC(1);
            Xrefined(2*i,3) = C(1);
            Yrefined(2*i-1,1) = A(2);
            Yrefined(2*i,1) = A(2);
            Yrefined(2*i-1,2) = B(2);
            Yrefined(2*i,2) = midptBC(2);
            Yrefined(2*i-1,3) = midptBC(2);
            Yrefined(2*i,3) = C(2);
            Zrefined(2*i-1,1) = A(3);
            Zrefined(2*i,1) = A(3);
            Zrefined(2*i-1,2) = B(3);
            Zrefined(2*i,2) = midptBC(3);
            Zrefined(2*i-1,3) = midptBC(3);
            Zrefined(2*i,3) = C(3);
        elseif maxInd ==3
            midptAC = A + 0.5*AC;
            Xrefined(2*i-1,1) = A(1);
            Xrefined(2*i,1) = midptAC(1);
            Xrefined(2*i-1,2) = B(1);
            Xrefined(2*i,2) = B(1);
            Xrefined(2*i-1,3) = midptAC(1);
            Xrefined(2*i,3) = C(1);
            Yrefined(2*i-1,1) = A(2);
            Yrefined(2*i,1) = midptAC(2);
            Yrefined(2*i-1,2) = B(2);
            Yrefined(2*i,2) = B(2);
            Yrefined(2*i-1,3) = midptAC(2);
            Yrefined(2*i,3) = C(2);
            Zrefined(2*i-1,1) = A(3);
            Zrefined(2*i,1) = midptAC(3);
            Zrefined(2*i-1,2) = B(3);
            Zrefined(2*i,2) = B(3);
            Zrefined(2*i-1,3) = midptAC(3);
            Zrefined(2*i,3) = C(3);
        end
    end
    X = Xrefined;
    Y = Yrefined;
    Z = Zrefined;
end
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

% Example of manual refinement
% this function splits the triangles in half N times, where N is the
% 4th argument of the function
% [X,Y,Z] = refineXYZ(X,Y,Z,1)
planes = [X(:,1) Y(:,1) Z(:,1) X(:,2) Y(:,2) Z(:,2) X(:,3) Y(:,3) Z(:,3)];
figure(3)
hold on
patch(transpose(X),transpose(Y),transpose(Z),color,'FaceAlpha',1,'EdgeColor','k')%'none')
end