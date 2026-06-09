close all
clear all

% iter_profile_solps;
% clear var X Y;
nP = 1e6;

% Load STL Geometry
model = createpde;
importGeometry(model,'west.stl');
figure(2)
pdegplot(model,'FaceLabels','on')

tic
mesh = generateMesh(model,'GeometricOrder','linear','Hmax',0.5);
figure(1)
pdeplot3D(model,'FaceAlpha',0.5)
[p,e,t] = meshToPet(mesh);

tess = transpose(t(1:4,:));
faces = [tess(:,[1 2 3]);tess(:,[1 2 4]); tess(:,[1 3 4]);tess(:,[2 3 4])];
faces = sort(faces,2);
faces = sortrows(faces);
Y = diff(faces);
zeroRow = [0,0,0];
k = ismember(Y,zeroRow,'rows');
k2 = find(k~=0);
faces([k2;k2+1],:) = [];

C = faces;
planes1 = [transpose(p(1:3,C(:,1))), transpose(p(1:3,C(:,2))), transpose(p(1:3,C(:,3)))];
% Adjust the coordinates a/to geometry
planes1_copy=planes1;
planes1(:,2) = planes1_copy(:,3);
planes1(:,5) = planes1_copy(:,6);
planes1(:,8) = planes1_copy(:,9);
planes1(:,3) = planes1_copy(:,2);
planes1(:,6) = planes1_copy(:,5);
planes1(:,9) = planes1_copy(:,8);


% Load second geometry from .mat file
load('west_antenna_geom.mat');  % Assumes 'planes' variable is present

% Combine STL + MAT geometry
planes = [planes1; planes];  % Replace 'planes' if needed
toc

% Identify which triangles came from final_ITERGeom
index_mat_start = size(planes1, 1) + 1;
index_mat_end = size(planes, 1);
inds_mat = index_mat_start:index_mat_end;

% Compute centroids
centroid = [1/3*(planes(:,1)+planes(:,4)+planes(:,7)), ...
            1/3*(planes(:,2)+planes(:,5)+planes(:,8)), ...
            1/3*(planes(:,3)+planes(:,6)+planes(:,9))];

nFaces = size(planes,1);
area = zeros(nFaces,1);
abcd = zeros(nFaces,4);
plane_norm = zeros(nFaces,1);
BCxBA = zeros(nFaces,1);
CAxCB = zeros(nFaces,1);

for i = 1:nFaces
    A = planes(i,1:3); B = planes(i,4:6); C = planes(i,7:9);
    AB = B - A; AC = C - A; BC = C - B; BA = A - B; CA = -AC; CB = -BC;

    norm1 = norm(AB); norm2 = norm(BC); norm3 = norm(AC);
    s = (norm1 + norm2 + norm3)/2;
    area(i) = sqrt(s*(s - norm1)*(s - norm2)*(s - norm3));
    normalVec = cross(AB, AC);
    abcd(i,:) = [normalVec, -dot(normalVec, A)];
    plane_norm(i) = norm(normalVec);
    BCxBA(i) = sign(dot(cross(BC, BA), normalVec));
    CAxCB(i) = sign(dot(cross(CA, CB), normalVec));
end

% Initialize potential and material fields
materialZ = 10*zeros(nFaces,1);
surfs = ones(nFaces,1);
potential_surf = zeros(nFaces, 1);
inDir = ones(1,length(planes));

%  dot_prod4selectSurf=0*inDir;
norm_vec=ones(length(surfs),3);
for i=1:length(surfs)
    ind = i;
    normal = -abcd(ind,1:3)/plane_norm(ind);
    l_normal = 0.1;
    normal = l_normal*normal;

    % dot_product = dot(normal,face_to_center(i,:));
    % 
    % if(sign(dot_product) == 1)
    %     inDir(ind) = -1;
    % end

    % dot_prod4selectSurf(ind)=dot(inDir(ind).*normal,[-2.6,1.4,0]);

    % 
    % if( dot_prod4selectSurf(ind) < -0.2)
    %     surfs(ind) = 1;
    %     materialZ(ind) = 74;
    % end


    normal = inDir(ind)*normal;
    norm_vec(i,:)=normal;
    vec_factor = 1;
    % quiver3(centroid(i,1),centroid(i,2),centroid(i,3),vec_factor*normal(1),vec_factor*normal(2),vec_factor*normal(3),'r')

end

% Write profiles for ITER

% save('final_ITERGeom.mat')


% Load field data
profilesWEST
antenna_data_WEST

%% GITR geometry

% Sheath calculations
disp('>>>> Calculating surface potentials')

comsol_data=readmatrix('../comsol/VDC_out_WEST57877_it3.txt');
xx = comsol_data(:,1);
yy = comsol_data(:,2);
zz = comsol_data(:,3);
emag = comsol_data(:,4);

% Apply potential only to triangles from .mat geometry
for i = inds_mat
     distance = sqrt((centroid(i,1) - xx).^2 + (centroid(i,2) - yy).^2 + (centroid(i,3) - zz).^2);

            [M I] = min(distance);
            potential_surf(i) = emag(I);
       
        end
         potential_surf=potential_surf';



% Plot for visualization
X = [planes(:,1), planes(:,4), planes(:,7)];
Y = [planes(:,2), planes(:,5), planes(:,8)];
Z = [planes(:,3), planes(:,6), planes(:,9)];
figure;
patch(transpose(X), transpose(Y), transpose(Z), potential_surf, 'FaceAlpha', 0.1, 'EdgeColor', 'none');
title('Potential applied only on ITER geometry')
xlabel('X'); ylabel('Y'); zlabel('Z');

figure; 
scatter3(centroid(:,1),centroid(:,2),centroid(:,3), 20, potential_surf, 'filled');
% Continue with particle generation, netCDF, etc., as in your full pipeline...

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


% Compute erosion flux only on .mat triangles
data = readmatrix('Targets.txt');
yields_data = readmatrix('yields_ITER_Ne2.csv');
% yields_data = 0.01;

ne_data = data(:,2);
v_data = data(:,5);
erosion_flux = zeros(nFaces, 1);
erosion_flux(inds_mat) = [0; yields_data(:,2)] .* ne_data .* v_data;
erosion_flux(isnan(erosion_flux)) = 0;
erosion = erosion_flux .* area;
% data=readmatrix('Targets.txt');
% yields_data=readmatrix('yields_ITER_Ne2.csv');
% 
% potential_data=data(:,1);
% ne_data=data(:,2);
% te_data=data(:,3);
% v_data=data(:,5);
% erosion_flux=[0;yields_data(:,2)].*ne_data.*v_data;
% erosion_flux(isnan(erosion_flux))=0;




% Erosion data
% Y0 = readmatrix('yields.csv');
% eroded_flux= ; % proportional to density now, need some realistic values
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
function ret = plot_panes(planes0,plotSet)

X = [planes0((plotSet),1),planes0((plotSet),4),planes0((plotSet),7)];
Y = [planes0((plotSet),2),planes0((plotSet),5),planes0((plotSet),8)];
Z = [planes0((plotSet),3),planes0((plotSet),6),planes0((plotSet),9)];

% [X,Y,Z] = refineXYZ(X,Y,Z,6)
planes0 = [X(:,1) Y(:,1) Z(:,1) X(:,2) Y(:,2) Z(:,2) X(:,3) Y(:,3) Z(:,3)];
figure
patch(transpose(X),transpose(Y),transpose(Z),[.7 .7 .7],'FaceAlpha',.3,'EdgeColor','k')%'none')


title({'GITR Geometry'})
xlabel('X [m]')
ylabel('Y [m]')
zlabel('Z [m]')
ret = 1;
end

function planes = planes_from_stl(filename,hmax)
model = createpde;
importGeometry(model,filename);% Import STL file
figure(2)
pdegplot(model,'FaceLabels','on') %Plot stl 

tic %time meshing - for high resolution this is a large cost
mesh = generateMesh(model,'GeometricOrder','linear','Hmax',hmax);% Options Hmax and Hmin can be set, linear order can also be used
figure
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

end

function planes_refined = refine_planes(planes,n)
X = [planes(:,1),planes(:,4),planes(:,7)];
Y = [planes(:,2),planes(:,5),planes(:,8)];
Z = [planes(:,3),planes(:,6),planes(:,9)];

[Xrefined, Yrefined, Zrefined] = refineXYZ(X,Y,Z,n);

planes_refined = [Xrefined(:,1) Yrefined(:,1) Zrefined(:,1) Xrefined(:,2) Yrefined(:,2) Zrefined(:,2) Xrefined(:,3) Yrefined(:,3) Zrefined(:,3)];
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
